from Bio import Entrez
import pandas as pd
import re
import os
import time
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlparse, parse_qs, unquote
from datetime import datetime

# --------------------------------------------------
# NCBI SETTINGS
# --------------------------------------------------
# Set NCBI_API_KEY as an environment variable for ~3x faster extraction
# (free key: https://account.ncbi.nlm.nih.gov/settings/ -> API Key Management)

Entrez.email = os.environ.get("NCBI_EMAIL", "bgpeeru@gmail.com")
API_KEY = os.environ.get("NCBI_API_KEY")
if API_KEY:
    Entrez.api_key = API_KEY

# NCBI limits: 3 requests/sec without a key, 10/sec with a key
if API_KEY:
    MAX_WORKERS = 8
    MIN_INTERVAL = 0.11
else:
    MAX_WORKERS = 4
    MIN_INTERVAL = 0.35

BATCH_SIZE = 200
MAX_RETRIES = 4
RETRY_DELAY = 1.5
MAX_RECORDS = 9999  # esearch hard limit per query

OUTPUT_DIR = os.path.join(tempfile.gettempdir(), "pubmed_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --------------------------------------------------
# RATE LIMITER (shared by all threads)
# --------------------------------------------------

_rate_lock = threading.Lock()
_last_call = [0.0]


def _throttle():
    with _rate_lock:
        wait = MIN_INTERVAL - (time.monotonic() - _last_call[0])
        if wait > 0:
            time.sleep(wait)
        _last_call[0] = time.monotonic()


# --------------------------------------------------
# COUNTRY DETECTION
# --------------------------------------------------

COUNTRY_LIST = [
    "United States of America", "United States", "USA", "U.S.A.",
    "United Kingdom", "UK", "England",
    "Scotland", "Wales", "Ireland", "Canada", "Australia",
    "People's Republic of China", "PR China", "China", "Japan", "India",
    "Germany", "France", "Italy", "Spain", "Netherlands",
    "Belgium", "Switzerland", "Austria", "Sweden", "Norway",
    "Denmark", "Finland", "Poland", "Brazil", "Mexico",
    "Singapore", "Republic of Korea", "South Korea", "Korea", "Taiwan",
    "Hong Kong", "Israel", "Turkey", "Portugal",
    "New Zealand", "Saudi Arabia", "UAE", "Qatar",
    "South Africa", "Russia", "Thailand", "Malaysia"
]

COUNTRY_MAP = {
    "USA": "United States",
    "U.S.A.": "United States",
    "United States of America": "United States",
    "UK": "United Kingdom",
    "England": "United Kingdom",
    "Scotland": "United Kingdom",
    "Wales": "United Kingdom",
    "People's Republic of China": "China",
    "PR China": "China",
    "Republic of Korea": "South Korea",
    "Korea": "South Korea"
}

# Single compiled regex with word boundaries (avoids "uk" matching inside
# words like "Luke", and is much faster than looping over the list)
_COUNTRY_RE = re.compile(
    r"(?<![A-Za-z])(" + "|".join(re.escape(c) for c in COUNTRY_LIST) + r")(?![A-Za-z])",
    re.IGNORECASE
)
_COUNTRY_LOOKUP = {c.lower(): COUNTRY_MAP.get(c, c) for c in COUNTRY_LIST}

_EMAIL_RE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")

_TRIAL_RE = re.compile(
    r"NCT\d{8}|ISRCTN\d+|ACTRN\d+|ChiCTR\d+|UMIN\d+|"
    r"CTRI/\d{4}/\d+/\d+|\d{4}-\d{6}-\d{2}",
    re.IGNORECASE
)


def extract_country(affiliation):
    if not affiliation:
        return ""
    m = _COUNTRY_RE.search(affiliation)
    if not m:
        return ""
    return _COUNTRY_LOOKUP.get(m.group(1).lower(), "")


def extract_trial_ids(article):
    return "; ".join(sorted(set(_TRIAL_RE.findall(str(article)))))


def safe_text(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return " ".join(str(x) for x in value)
    return str(value)


# --------------------------------------------------
# NCBI CALLS WITH RETRY
# --------------------------------------------------

def _with_retry(func, label, **kwargs):
    last_error = None
    for attempt in range(1, MAX_RETRIES + 1):
        handle = None
        try:
            _throttle()
            handle = func(**kwargs)
            return Entrez.read(handle)
        except Exception as e:
            last_error = e
            print(f"{label} failed (attempt {attempt}/{MAX_RETRIES}): {e}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY * attempt)
        finally:
            if handle is not None:
                try:
                    handle.close()
                except Exception:
                    pass
    raise RuntimeError(f"{label} failed after {MAX_RETRIES} attempts: {last_error}")


def search_pubmed_with_retry(query, retmax):
    return _with_retry(
        Entrez.esearch, "Search",
        db="pubmed", term=query, retmax=retmax, retmode="xml"
    )


def fetch_batch(batch_ids):
    return _with_retry(
        Entrez.efetch, "Fetch",
        db="pubmed", id=",".join(batch_ids), retmode="xml"
    )


# --------------------------------------------------
# PARSE ONE ARTICLE -> LIST OF AUTHOR ROWS
# --------------------------------------------------

def parse_article(article):
    rows = []

    medline = article["MedlineCitation"]
    article_data = medline["Article"]
    pmid = str(medline["PMID"])

    title = safe_text(article_data.get("ArticleTitle", ""))

    journal = ""
    year = ""
    try:
        journal = safe_text(article_data["Journal"]["Title"])
        pub_date = article_data["Journal"]["JournalIssue"]["PubDate"]
        year = safe_text(pub_date.get("Year", "")) or safe_text(
            pub_date.get("MedlineDate", "")
        )
    except Exception:
        pass

    doi = ""
    try:
        for eid in article_data.get("ELocationID", []):
            if eid.attributes.get("EIdType") == "doi":
                doi = str(eid)
                break
    except Exception:
        pass

    abstract = ""
    try:
        parts = []
        for item in article_data["Abstract"]["AbstractText"]:
            label = item.attributes.get("Label", "")
            parts.append(f"{label}: {item}" if label else str(item))
        abstract = " ".join(parts)
    except Exception:
        pass

    trial_id = extract_trial_ids(article)

    authors = article_data.get("AuthorList", [])
    total_authors = len(authors)

    for i, author in enumerate(authors):
        try:
            if "ForeName" not in author:
                continue

            fullname = (
                safe_text(author.get("ForeName", ""))
                + " "
                + safe_text(author.get("LastName", ""))
            ).strip()

            if i == 0:
                position = "First"
            elif i == total_authors - 1:
                position = "Last"
            else:
                position = "Middle"

            affiliation = ""
            email = ""
            country = ""

            if "AffiliationInfo" in author:
                aff_list = []
                for aff in author["AffiliationInfo"]:
                    txt = safe_text(aff.get("Affiliation", ""))
                    aff_list.append(txt)
                    if not email:
                        m = _EMAIL_RE.search(txt)
                        if m:
                            email = m.group()
                affiliation = "; ".join(aff_list)
                country = extract_country(affiliation)

            rows.append({
                "PMID": pmid,
                "Title": title,
                "Journal": journal,
                "Year": year,
                "DOI": doi,
                "Author Full Name": fullname,
                "Author Position": position,
                "Affiliation": affiliation,
                "Country": country,
                "Email ID": email,
                "Trial ID": trial_id,
                "Grant ID": "",
                "Publication Type": "",
                "MeSH Terms": "",
                "Keywords": "",
                "Abstract": abstract,
                "Department": "",
                "Institution": ""
            })
        except Exception as e:
            print("Skipped author:", e)

    return rows


def process_batch(batch_ids):
    records = fetch_batch(batch_ids)
    rows = []
    for article in records.get("PubmedArticle", []):
        try:
            rows.extend(parse_article(article))
        except Exception as e:
            print("Skipped article:", e)
    return rows


# --------------------------------------------------
# MAIN EXTRACTOR
# --------------------------------------------------

def extract_pubmed(
        search_type,
        keyword=None,
        start_year=None,
        end_year=None,
        url=None,
        retmax=500,
        output_choice="1",
        selected_columns=None
):
    t0 = time.time()

    if selected_columns is None:
        selected_columns = []

    retmax = max(1, min(int(retmax), MAX_RECORDS))

    # ---------------- BUILD QUERY ----------------
    if search_type == "keyword":
        QUERY = (
            f'("{keyword}"[Title/Abstract]) '
            f'AND ("{start_year}/01/01"[Date - Publication] : '
            f'"{end_year}/12/31"[Date - Publication])'
        )
    elif search_type == "url":
        parsed = urlparse(url)
        params = parse_qs(parsed.query)
        if "term" not in params:
            raise ValueError("Invalid PubMed URL")
        QUERY = unquote(params["term"][0])
        if "filter" in params:
            for f in params["filter"]:
                if f == "articleattr.data":
                    QUERY += " AND data[sb]"
    else:
        raise ValueError("Invalid Search Type")

    print("Query:", QUERY, "| Requested:", retmax)

    # ---------------- SEARCH ----------------
    search = search_pubmed_with_retry(QUERY, retmax)
    pmids = list(search["IdList"])
    print("Total found:", search.get("Count", 0), "| Selected:", len(pmids))

    if not pmids:
        raise Exception("No records found.")

    # ---------------- PARALLEL FETCH ----------------
    batches = [pmids[i:i + BATCH_SIZE] for i in range(0, len(pmids), BATCH_SIZE)]
    results = {}
    failed = []

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(process_batch, b): idx for idx, b in enumerate(batches)}
        for fut in as_completed(futures):
            idx = futures[fut]
            try:
                results[idx] = fut.result()
                print(f"Batch {idx + 1}/{len(batches)} done ({len(results[idx])} rows)")
            except Exception as e:
                failed.append(idx)
                print(f"Batch {idx + 1} FAILED: {e}")

    # Keep original PubMed order
    rows = []
    for idx in sorted(results):
        rows.extend(results[idx])

    df = pd.DataFrame(rows)

    if df.empty:
        raise Exception(
            "No records could be extracted. "
            "Check the search query or PubMed response."
        )

    df = df.drop_duplicates(subset=["PMID", "Author Full Name"])

    # ---------------- OUTPUT SELECTION ----------------
    is_custom_output = str(output_choice).strip().lower() in {
        "2", "custom", "custom columns"
    }

    if is_custom_output:
        if not selected_columns:
            raise Exception("Please select at least one column.")

        selected_columns = list(dict.fromkeys(selected_columns))
        missing = [c for c in selected_columns if c not in df.columns]
        if missing:
            raise Exception(f"Invalid columns selected: {missing}")

        final_df = df.loc[:, selected_columns].drop_duplicates()
    else:
        final_df = df

    # ---------------- SAVE ----------------
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(OUTPUT_DIR, f"PubMed_Author_Level_{timestamp}.csv")
    final_df.to_csv(output_file, index=False, encoding="utf-8-sig")

    print(
        f"DONE: {len(final_df)} rows, {len(batches) - len(failed)}/{len(batches)} "
        f"batches OK, {time.time() - t0:.1f}s -> {output_file}"
    )

    return output_file