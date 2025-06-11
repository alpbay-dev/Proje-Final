# app/routes.py

from flask import Blueprint, render_template, request, url_for
import os
import shutil

from app.helpers.Duru.boldsystems_api import fetch_sequences_by_organism_and_marker
from app.helpers.Duru.core.sequence_processing import create_fasta_from_sequences, SequenceProcessingError
from app.helpers.Tuana.AnalysisFunctions import (
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

@main.route("/aboutus")
def aboutus():
    return render_template("aboutus.html")


@main.route("/align", methods=["GET", "POST"])
def align():
    if request.method == "POST":
        print(">>> POST isteği alındı.")
        marker_code = request.form["dna_type"]
        species_list = [s.strip() for s in request.form.getlist("species[]") if s.strip()]
        seq_limit = int(request.form["sequence_count"])
        selected_analyses = request.form.getlist("analysis")

        print(f">>> Marker: {marker_code}")
        print(f">>> Tür listesi: {species_list}")
        print(f">>> Sekans limiti: {seq_limit}")
        print(f">>> Seçilen analizler: {selected_analyses}")

        all_sequences = []
        warnings = []

        for sp in species_list:
            try:
                print(f">>> BOLD API çağrılıyor: {sp}")
                seqs = fetch_sequences_by_organism_and_marker(
                    organism=sp,
                    marker_code=marker_code,
                    limit=seq_limit
                )
                if len(seqs) < seq_limit:
                    warnings.append(f"{sp}: Requested {seq_limit} sequences, but only {len(seqs)} were found.")
                all_sequences.extend(seqs)
            except Exception as e:
                print(f"❌ {sp} için veri çekme hatası: {e}")
                warnings.append(f"{sp}: Failed to fetch sequences. ({str(e)})")

        print(f">>> Toplam çekilen sekans sayısı: {len(all_sequences)}")

        if not all_sequences:
            return render_template("align.html", error="No sequences found for the selected species and marker.")

        os.makedirs(TEMP_FASTA_DIR, exist_ok=True)
        temp_fasta = os.path.join(TEMP_FASTA_DIR, "user_sequences.fasta")

        try:
            print(">>> FASTA oluşturuluyor")
            final_fasta = create_fasta_from_sequences(all_sequences, temp_fasta, align=True)
            print(">>> Aligned FASTA oluşturuldu")
        except SequenceProcessingError as e:
            print(f"❌ FASTA oluşturma hatası: {e}")
            return render_template("align.html", error=f"Failed to create FASTA: {e}")

        results = {}

        def save_result(title, image_path):
            if not os.path.exists(image_path):
                warnings.append(f"Unable to display result image for: {title}")
                return
            filename = os.path.basename(image_path)
            static_path = os.path.join(STATIC_IMG_DIR, filename)
            shutil.copyfile(image_path, static_path)
            results[title] = url_for("static", filename=f"analysis_results/{filename}")

        # Analizler
        if "ml" in selected_analyses:
            try:
                ml_path = run_and_draw_bootstrap_tree(
                    final_fasta,
                    iqtree_path=r"C:\Users\SARP\Desktop\Alp\Proje-Final\app\helpers\Tuana\iqtree\iqtree3.exe"
                )
                save_result("Maximum Likelihood (IQ-TREE)", ml_path)
            except Exception as e:
                warnings.append(f"Maximum Likelihood Tree failed: {e}")

        # Aligned FASTA her zaman eklenecek
        try:
            aligned_fasta_filename = os.path.basename(final_fasta)
            static_aligned_path = os.path.join(STATIC_IMG_DIR, aligned_fasta_filename)
            shutil.copyfile(final_fasta, static_aligned_path)
            results["Aligned FASTA File"] = url_for("static", filename=f"analysis_results/{aligned_fasta_filename}")
        except Exception as e:
            warnings.append("Aligned FASTA file could not be found.")

        if "msn" in selected_analyses:
            try:
                msn_path = build_msn_with_distancecalculator(final_fasta)
                save_result("Minimum Spanning Network", msn_path)
            except Exception as e:
                warnings.append(f"Minimum Spanning Network failed: {e}")

        if "heatmap" in selected_analyses:
            try:
                heatmap_path = os.path.splitext(final_fasta)[0] + "_pi_heatmap.png"
                plot_nucleotide_diversity_heatmap(final_fasta, output_img=heatmap_path)
                save_result("Nucleotide Diversity Heatmap", heatmap_path)
            except Exception as e:
                warnings.append(f"Nucleotide Diversity Heatmap failed: {e}")

        return render_template("results.html", result=results, warnings=warnings)

    return render_template("align.html", result=None)
