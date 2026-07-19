from Bio import Entrez
import pandas as pd
from tqdm import tqdm
from urllib.parse import urlparse, parse_qs, unquote
import re
import time
import os

Entrez.email = "bgpeeru@gmail.com"

BATCH_SIZE = 100


def extract_pubmed(search_type,
                   keyword=None,
                   start_year=None,
                   end_year=None,
                   url=None,
                   retmax=500):
    """
    Extract PubMed author details and save them to CSV.
    Returns the path of the generated CSV.
    """

    # -----------------------------
    # Build PubMed Query
    # -----------------------------
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

    else:
        raise ValueError("Invalid Search Type")

    print("Query:", QUERY)

    # -----------------------------
    # Search PubMed
    # -----------------------------
    handle = Entrez.esearch(
        db="pubmed",
        term=QUERY,
        retmax=retmax
    )

    search = Entrez.read(handle)
    handle.close()

    pmids = search["IdList"]

    print("Total PMIDs:", len(pmids))

    rows = []

    # -----------------------------
    # Fetch Records
    # -----------------------------
    for start in tqdm(range(0, len(pmids), BATCH_SIZE)):

        end = min(start + BATCH_SIZE, len(pmids))

        batch_ids = ",".join(pmids[start:end])

        handle = Entrez.efetch(
            db="pubmed",
            id=batch_ids,
            retmode="xml"
        )

        records = Entrez.read(handle)
        handle.close()

        for article in records["PubmedArticle"]:

            try:

                medline = article["MedlineCitation"]
                article_data = medline["Article"]

                pmid = str(medline["PMID"])

                year = ""

                try:
                    year = article_data["Journal"]["JournalIssue"]["PubDate"]["Year"]
                except:
                    pass

                authors = article_data.get("AuthorList", [])

                total = len(authors)

                for i, author in enumerate(authors):

                    if "ForeName" not in author:
                        continue

                    fullname = (
                        author.get("ForeName", "")
                        + " "
                        + author.get("LastName", "")
                    ).strip()

                    if i == 0:
                        position = "First"
                    elif i == total - 1:
                        position = "Last"
                    else:
                        position = "Middle"

                    affiliation = ""
                    email = ""

                    if "AffiliationInfo" in author:

                        affs = []

                        for aff in author["AffiliationInfo"]:

                            txt = aff.get("Affiliation", "")

                            affs.append(txt)

                            m = re.search(
                                r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}",
                                txt
                            )

                            if m:
                                email = m.group()

                        affiliation = "; ".join(affs)

                    rows.append({
                        "PMID": pmid,
                        "Author Full Name": fullname,
                        "Author Position": position,
                        "Affiliation": affiliation,
                        "Email ID": email,
                        "Year": year
                    })

            except Exception:
                continue

        time.sleep(0.34)

    # -----------------------------
    # Create DataFrame
    # -----------------------------
    df = pd.DataFrame(rows)

    df = df.drop_duplicates(
        subset=["PMID", "Author Full Name"]
    )

    # -----------------------------
    # Create Output Folder
    # -----------------------------
    os.makedirs("output", exist_ok=True)

    output_file = os.path.join(
        "output",
        "PubMed_Author_Level.csv"
    )

    df.to_csv(output_file, index=False)

    print("Completed Successfully!")
    print("Author Records:", len(df))

    return output_file