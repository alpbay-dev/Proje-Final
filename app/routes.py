from flask import Blueprint, render_template, request
from app.helpers.one import process_data_one
from app.helpers.two import process_data_two

main = Blueprint("main", __name__)

@main.route("/")
def home():
    return render_template("home.html")

@main.route("/align", methods=["GET", "POST"])
def align():
    result = None
    if request.method == "POST":
        dna_type = request.form.get("dna_type")
        species = request.form.getlist("species[]")
        seq_count = request.form.get("seq_count")
        geo_div = request.form.get("geo_div") == "on"

        # Şimdilik sadece işlemleri yazdırıyoruz
        result_one = process_data_one(dna_type, species)
        result_two = process_data_two(seq_count, geo_div)

        result = f"ONE: {result_one}<br>TWO: {result_two}"

    return render_template("align.html", result=result)
