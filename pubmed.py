from Bio import Entrez
import pandas as pd
import re
import time
from tqdm import tqdm
from urllib.parse import urlparse, parse_qs, unquote

Entrez.email = "your_email@gmail.com"

BATCH_SIZE = 100

COUNTRY_LIST = [
    "United States","USA","United Kingdom","UK","England",
    "Scotland","Wales","Ireland","Canada","Australia",
    "China","People's Republic of China","Japan","India",
    "Germany","France","Italy","Spain","Netherlands",
    "Belgium","Switzerland","Austria","Sweden","Norway",
    "Denmark","Finland","Poland","Brazil","Mexico",
    "Singapore","South Korea","Korea","Taiwan",
    "Hong Kong","Israel","Turkey","Portugal",
    "New Zealand","Saudi Arabia","UAE","Qatar",
    "South Africa","Russia","Thailand","Malaysia"
]

COUNTRY_MAP = {
    "USA":"United States",
    "U.S.A.":"United States",
    "United States of America":"United States",
    "UK":"United Kingdom",
    "England":"United Kingdom",
    "Scotland":"United Kingdom",
    "Wales":"United Kingdom",
    "People's Republic of China":"China",
    "PR China":"China",
    "Republic of Korea":"South Korea",
    "Korea":"South Korea"
}


def extract_country(affiliation):

    if not affiliation:
        return ""

    aff = affiliation.lower()

    for country in COUNTRY_LIST:

        if country.lower() in aff:
            return COUNTRY_MAP.get(country, country)

    return ""


def extract_trial_ids(article):

    trial_ids = []

    text = str(article)

    patterns = [
        r"NCT\d{8}",
        r"ISRCTN\d+",
        r"ACTRN\d+",
        r"ChiCTR\d+",
        r"UMIN\d+",
        r"CTRI/\d{4}/\d+/\d+",
        r"\d{4}-\d{6}-\d{2}"
    ]

    for pattern in patterns:

        ids = re.findall(pattern, text)

        if ids:
            trial_ids.extend(ids)

    return "; ".join(sorted(set(trial_ids)))


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

    if selected_columns is None:
        selected_columns = []

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

    handle = Entrez.esearch(
        db="pubmed",
        term=QUERY,
        retmax=retmax
    )

    search = Entrez.read(handle)

    handle.close()

    pmids = search["IdList"]

    rows = []
    for start in tqdm(range(0, len(pmids), BATCH_SIZE), desc="Processing Articles"):
        end = min(start + BATCH_SIZE, len(pmids))
        batch_ids = ",".join(pmids[start:end])

        handle = Entrez.efetch(db="pubmed", id=batch_ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for article in records["PubmedArticle"]:
            try:
                medline = article["MedlineCitation"]
                article_data = medline["Article"]
                pmid = str(medline["PMID"])

                title = ""
                try:
                    title = article_data["ArticleTitle"]
                except Exception:
                    pass

                journal = ""
                try:
                    journal = article_data["Journal"]["Title"]
                except Exception:
                    pass

                year = ""
                try:
                    year = article_data["Journal"]["JournalIssue"]["PubDate"]["Year"]
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
                    abstract = " ".join(article_data["Abstract"]["AbstractText"])
                except Exception:
                    pass

                trial_id = extract_trial_ids(article)
                authors = article_data.get("AuthorList", [])
                total = len(authors)

                for i, author in enumerate(authors):
                    if "ForeName" not in author:
                        continue

                    fullname = (
                        author.get("ForeName", "") + " " + author.get("LastName", "")
                    ).strip()

                    if i == 0:
                        position = "First"
                    elif i == total - 1:
                        position = "Last"
                    else:
                        position = "Middle"

                    affiliation = ""
                    email = ""
                    country = ""

                    if "AffiliationInfo" in author:
                        aff_list = []
                        for aff in author["AffiliationInfo"]:
                            txt = aff.get("Affiliation", "")
                            aff_list.append(txt)

                            match = re.search(
                                r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}",
                                txt,
                            )

                            if match:
                                email = match.group()

                        affiliation = "; ".join(aff_list)
                        country = extract_country(affiliation)

                    rows.append(
                        {
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
                            "Institution": "",
                        }
                    )
            except Exception as e:
                print("Skipped:", e)
                continue

        time.sleep(0.34)

    # ---------------------------------
    # Create DataFrame
    # ---------------------------------

    df = pd.DataFrame(rows)

    if df.empty:
        raise Exception("No records found.")

    # Remove duplicate author rows
    df = df.drop_duplicates(subset=["PMID", "Author Full Name"])

    # ---------------------------------
    # Custom Output
    # ---------------------------------

    if output_choice == "2":

        if not selected_columns:
            raise Exception("Please select at least one column.")

        missing = [c for c in selected_columns if c not in df.columns]

        if missing:
            raise Exception(f"Invalid columns selected: {missing}")

        final_df = df[selected_columns]

        final_df = final_df.drop_duplicates()

    else:

        final_df = df

    # ---------------------------------
    # Save CSV
    # ---------------------------------

    output_file = "PubMed_Author_Level.csv"

    final_df.to_csv(output_file, index=False)

    return output_file