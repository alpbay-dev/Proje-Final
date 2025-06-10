# app/routes.py

from flask import Blueprint, render_template, request, url_for
import os
import shutil


from app.helpers.Duru.boldsystems_api import fetch_sequences_by_organism_and_marker


from app.helpers.Duru.core.sequence_processing import create_fasta_from_sequences, SequenceProcessingError


from app.helpers.Tuana.AnalysisFunctions import (
    build_nj_tree,
    build_mp_tree,
    run_and_draw_bootstrap_tree,
    build_msn_with_distancecalculator,
    plot_nucleotide_diversity_heatmap
)

main = Blueprint("main", __name__)


BASE_DIR = os.path.dirname(__file__)
TEMP_FASTA_DIR = os.path.join(BASE_DIR, "helpers", "Duru", "data", "temp_fasta")
STATIC_IMG_DIR = os.path.join("app", "static", "analysis_results")
os.makedirs(TEMP_FASTA_DIR, exist_ok=True)
os.makedirs(STATIC_IMG_DIR, exist_ok=True)

@main.route("/")
def home():
    return render_template("home.html")


@main.route("/align", methods=["GET", "POST"])
def align():
    if request.method == "POST":

        marker_code = request.form["dna_type"]
        species_list = [s.strip() for s in request.form.getlist("species[]") if s.strip()]
        seq_limit = int(request.form["sequence_count"])
        selected_analyses = request.form.getlist("analysis")


        all_sequences = []
        for sp in species_list:
            try:
                seqs = fetch_sequences_by_organism_and_marker(
                    organism=sp,
                    marker_code=marker_code,
                    limit=seq_limit
                )
                all_sequences.extend(seqs)
            except Exception as e:
                print(f"Error fetching {sp}: {e}")

        if not all_sequences:
            return render_template("align.html", error="No sequences found for the selected species and marker.")


        os.makedirs(TEMP_FASTA_DIR, exist_ok=True)
        temp_fasta = os.path.join(TEMP_FASTA_DIR, "user_sequences.fasta")
        try:
            final_fasta = create_fasta_from_sequences(all_sequences, temp_fasta, align=True)  # 👈 DİKKAT: align=True
        except SequenceProcessingError as e:
            return render_template("align.html", error=f"Failed to create FASTA: {e}")


        results = {}

        def save_result(title, image_path):
            filename = os.path.basename(image_path)
            static_path = os.path.join(STATIC_IMG_DIR, filename)
            shutil.copyfile(image_path, static_path)
            results[title] = url_for('static', filename=f'analysis_results/{filename}')


        if "nj" in selected_analyses:
            nj_path = build_nj_tree(final_fasta)
            save_result("Neighbor Joining Tree", nj_path)

        if "mp" in selected_analyses:
            mp_path = build_mp_tree(final_fasta)
            save_result("Maximum Parsimony Tree", mp_path)

        if "ml" in selected_analyses:
            ml_path = run_and_draw_bootstrap_tree(final_fasta)
            save_result("Maximum Likelihood (IQ-TREE)", ml_path)

        if "msn" in selected_analyses:
            msn_path = build_msn_with_distancecalculator(final_fasta)
            save_result("Minimum Spanning Network", msn_path)

        if "heatmap" in selected_analyses:
            heatmap_path = plot_nucleotide_diversity_heatmap(final_fasta)
            save_result("Nucleotide Diversity Heatmap", heatmap_path)


        return render_template("results.html", result=results)


    return render_template("align.html", result=None)
