from Bio import Entrez
import csv
import gc
import re
import os
import time
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlparse, parse_qs, unquote
from datetime import datetime

# --------------------------------------------------
# NCBI SETTINGS
# --------------------------------------------------
# Set NCBI_API_KEY as an environment variable for faster extraction
# (free key: https://account.ncbi.nlm.nih.gov/settings/ -> API Key Management)

Entrez.email = os.environ.get("NCBI_EMAIL", "bgpeeru@gmail.com")
API_KEY = os.environ.get("NCBI_API_KEY")
if API_KEY:
    Entrez.api_key = API_KEY

# NCBI limits: 3 requests/sec without a key, 10/sec with a key.
# Worker count and batch size are kept small on purpose: every worker holds
# one parsed XML batch in memory, and Render's free instance has only 512 MB.
if API_KEY:
    MAX_WORKERS = 4
    MIN_INTERVAL = 0.11
else:
    MAX_WORKERS = 3
    MIN_INTERVAL = 0.35

BATCH_SIZE = 100
MAX_RETRIES = 4
RETRY_DELAY = 1.5
MAX_RECORDS = 9999  # esearch hard limit per query

OUTPUT_DIR = os.path.join(tempfile.gettempdir(), "pubmed_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

ALL_COLUMNS = [
    "PMID", "Title", "Journal", "Year", "DOI", "Author Full Name",
    "Author Position", "Affiliation", "Country", "Email ID", "Trial ID",
    "Grant ID", "Publication Type", "MeSH Terms", "Keywords", "Abstract",
    "Department", "Institution"
]

# Only one extraction runs at a time, so two users clicking at once
# cannot double the memory use.
_EXTRACT_SLOT = threading.Semaphore(1)

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
# PUBMED URL FILTERS -> SEARCH QUERY
# --------------------------------------------------

_PUBTYPE_MAP = {
    "clinicaltrial": "Clinical Trial",
    "randomizedcontrolledtrial": "Randomized Controlled Trial",
    "review": "Review",
    "systematicreview": "Systematic Review",
    "metaanalysis": "Meta-Analysis",
    "meta-analysis": "Meta-Analysis",
    "casereports": "Case Reports",
    "observationalstudy": "Observational Study",
    "guideline": "Guideline",
    "editorial": "Editorial",
    "letter": "Letter",
}


def url_filters_to_query(filters):
    """
    Convert the sidebar filters in a PubMed URL (&filter=...) into
    search clauses. Filters of the same kind (e.g. several article
    types) are OR-ed; different kinds are AND-ed.
    """
    pubtypes = []
    clauses = []

    for f in filters:
        f = unquote(f).strip()
        key, _, value = f.partition(".")
        key = key.lower()

        if f == "articleattr.data":
            clauses.append("data[sb]")

        elif key == "years" and "-" in value:
            a, b = value.split("-", 1)
            clauses.append(f'("{a}/01/01"[dp] : "{b}/12/31"[dp])')

        elif key == "dates" and "-" in value:
            a, b = value.split("-", 1)
            clauses.append(f'("{a}"[dp] : "{b}"[dp])')

        elif key == "pubt" and value.lower() in _PUBTYPE_MAP:
            pubtypes.append(f'"{_PUBTYPE_MAP[value.lower()]}"[pt]')

        elif key == "lang":
            clauses.append(f"{value}[lang]")

        elif key == "hum_ani" and value.lower() in {"humans", "animals"}:
            clauses.append(f"{value.lower()}[mh]")

        elif key == "sex" and value.lower() in {"male", "female"}:
            clauses.append(f"{value.lower()}[mh]")

        elif f == "hasabstract":
            clauses.append("hasabstract")

        elif f == "simsearch1.fha":
            clauses.append("hasabstract")

        elif f == "simsearch2.ffrft":
            clauses.append('"free full text"[sb]')

        elif f == "simsearch3.fft":
            clauses.append('"full text"[sb]')

        else:
            print("Unsupported URL filter ignored:", f)

    if pubtypes:
        clauses.append("(" + " OR ".join(pubtypes) + ")")

    return clauses


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
# `need` is the set of columns the user asked for (None = all columns).
# Expensive fields (abstract, trial IDs, affiliations) are skipped when
# they are not needed, which saves both time and memory.

def parse_article(article, need=None):
    def want(*cols):
        return need is None or any(c in need for c in cols)

    rows = []

    medline = article["MedlineCitation"]
    article_data = medline["Article"]
    pmid = str(medline["PMID"])

    title = safe_text(article_data.get("ArticleTitle", "")) if want("Title") else ""

    journal = ""
    year = ""
    if want("Journal", "Year"):
        try:
            journal = safe_text(article_data["Journal"]["Title"])
            pub_date = article_data["Journal"]["JournalIssue"]["PubDate"]
            year = safe_text(pub_date.get("Year", "")) or safe_text(
                pub_date.get("MedlineDate", "")
            )
        except Exception:
            pass

    doi = ""
    if want("DOI"):
        try:
            for eid in article_data.get("ELocationID", []):
                if eid.attributes.get("EIdType") == "doi":
                    doi = str(eid)
                    break
        except Exception:
            pass

    abstract = ""
    if want("Abstract"):
        try:
            parts = []
            for item in article_data["Abstract"]["AbstractText"]:
                label = item.attributes.get("Label", "")
                parts.append(f"{label}: {item}" if label else str(item))
            abstract = " ".join(parts)
        except Exception:
            pass

    trial_id = extract_trial_ids(article) if want("Trial ID") else ""

    need_aff = want("Affiliation", "Country", "Email ID")

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

            if need_aff and "AffiliationInfo" in author:
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


def process_batch(batch_ids, need=None):
    records = fetch_batch(batch_ids)
    rows = []
    for article in records.get("PubmedArticle", []):
        try:
            rows.extend(parse_article(article, need))
        except Exception as e:
            print("Skipped article:", e)
    del records
    return rows


def _cleanup_old_outputs(max_age_seconds=3600):
    """Delete CSVs from earlier runs so the temp disk does not fill up."""
    now = time.time()
    try:
        for name in os.listdir(OUTPUT_DIR):
            path = os.path.join(OUTPUT_DIR, name)
            if os.path.isfile(path) and now - os.path.getmtime(path) > max_age_seconds:
                try:
                    os.remove(path)
                except OSError:
                    pass
    except OSError:
        pass


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

    # ---------------- OUTPUT SELECTION (validated up front) ----------------
    is_custom_output = str(output_choice).strip().lower() in {
        "2", "custom", "custom columns"
    }

    if is_custom_output:
        if not selected_columns:
            raise Exception("Please select at least one column.")

        columns = list(dict.fromkeys(selected_columns))
        missing = [c for c in columns if c not in ALL_COLUMNS]
        if missing:
            raise Exception(f"Invalid columns selected: {missing}")

        need = set(columns)
    else:
        columns = list(ALL_COLUMNS)
        need = None

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
        for clause in url_filters_to_query(params.get("filter", [])):
            QUERY += f" AND {clause}"
    else:
        raise ValueError("Invalid Search Type")

    print("Query:", QUERY, "| Requested:", retmax)

    # Only one extraction at a time (keeps memory use predictable)
    with _EXTRACT_SLOT:

        _cleanup_old_outputs()

        # ---------------- SEARCH ----------------
        search = search_pubmed_with_retry(QUERY, retmax)
        pmids = list(search["IdList"])
        print("Total found:", search.get("Count", 0), "| Selected:", len(pmids))

        if not pmids:
            raise Exception("No records found.")

        batches = [pmids[i:i + BATCH_SIZE] for i in range(0, len(pmids), BATCH_SIZE)]

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = os.path.join(OUTPUT_DIR, f"PubMed_Author_Level_{timestamp}.csv")

        def run(batch):
            try:
                return process_batch(batch, need)
            except Exception as e:
                print("Batch FAILED:", e)
                return None

        failed = 0
        total_rows = 0
        seen_authors = set()   # (PMID, author) pairs already written
        seen_rows = set()      # hashes of selected-column rows (custom output)

        # ---------------- PARALLEL FETCH, STREAMED TO DISK ----------------
        # Rows are written to the CSV as each batch finishes instead of being
        # collected in memory first. pool.map returns results in PubMed order.
        with open(output_file, "w", newline="", encoding="utf-8-sig") as f, \
                ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:

            writer = csv.writer(f)
            writer.writerow(columns)

            for idx, rows in enumerate(pool.map(run, batches), start=1):
                if rows is None:
                    failed += 1
                    continue

                for r in rows:
                    key = (r["PMID"], r["Author Full Name"])
                    if key in seen_authors:
                        continue
                    seen_authors.add(key)

                    values = [r[c] for c in columns]

                    if is_custom_output:
                        h = hash(tuple(values))
                        if h in seen_rows:
                            continue
                        seen_rows.add(h)

                    writer.writerow(values)
                    total_rows += 1

                print(f"Batch {idx}/{len(batches)} done | rows so far: {total_rows}")
                del rows
                gc.collect()

        if total_rows == 0:
            try:
                os.remove(output_file)
            except OSError:
                pass
            raise Exception(
                "No records could be extracted. "
                "Check the search query or PubMed response."
            )

    print(
        f"DONE: {total_rows} rows, {len(batches) - failed}/{len(batches)} "
        f"batches OK, {time.time() - t0:.1f}s -> {output_file}"
    )

    return output_file