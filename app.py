from flask import Flask, render_template, request, send_file
from pubmed import extract_pubmed

app = Flask(__name__)


@app.route("/")
def home():
    return render_template("index.html")


@app.route("/extract", methods=["POST"])
def extract():

    search_type = request.form.get("search_type")

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
                retmax=retmax
            )

        else:

            url = request.form.get("url")
            retmax = int(request.form.get("retmax"))

            csv_file = extract_pubmed(
                search_type="url",
                url=url,
                retmax=retmax
            )

        return send_file(csv_file, as_attachment=True)

    except Exception as e:

        return f"<h2>Error</h2><br>{str(e)}"


if __name__ == "__main__":
    app.run(debug=True)