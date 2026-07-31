from flask import Flask, render_template, request, send_file
from pubmed import extract_pubmed

app = Flask(__name__)


@app.route("/")
def home():
    return render_template("index.html")


@app.route("/extract", methods=["POST"])
def extract():

    search_type = request.form.get("search_type")
    output_choice = request.form.get("output_choice")
    selected_columns = request.form.getlist("columns")

    print("=" * 50)
    print("Search Type:", search_type)
    print("Output Choice:", output_choice)
    print("Selected Columns:", selected_columns)
    print("=" * 50)

    try:

        if search_type == "keyword":

            keyword = request.form.get("keyword")
            start_year = request.form.get("start_year")
            end_year = request.form.get("end_year")
            retmax = int(request.form.get("retmax"))

            csv_file = extract_pubmed(
                search_type="keyword",
                keyword=keyword,
                start_year=start_year,
                end_year=end_year,
                retmax=retmax,
                output_choice=output_choice,
                selected_columns=selected_columns
            )

        else:

            url = request.form.get("url")
            retmax = int(request.form.get("retmax"))

            csv_file = extract_pubmed(
                search_type="url",
                url=url,
                retmax=retmax,
                output_choice=output_choice,
                selected_columns=selected_columns
            )

        return send_file(csv_file, as_attachment=True)

    except Exception as e:

        return f"<h2>Error</h2><br>{str(e)}"


import os

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)