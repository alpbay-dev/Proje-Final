# app/routes.py

from flask import Blueprint, render_template, request, url_for
import os
import shutil

# 1. BOLD API
from app.helpers.Duru.boldsystems_api import fetch_sequences_by_organism_and_marker

# 2. FASTA işleme
from app.helpers.Duru.core.sequence_processing import create_fasta_from_sequences, SequenceProcessingError

# 3. Analiz fonksiyonları
from app.helpers.Tuana.AnalysisFunctions import (
    build_nj_tree,
    build_mp_tree,
    run_and_draw_bootstrap_tree,
    build_msn_with_distancecalculator,
    plot_nucleotide_diversity_heatmap
)

main = Blueprint("main", __name__)

# Klasör tanımlamaları
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
        print(">>> POST isteği alındı.")

        # Form verilerini al
        marker_code = request.form["dna_type"]
        species_list = [s.strip() for s in request.form.getlist("species[]") if s.strip()]
        seq_limit = int(request.form["sequence_count"])
        selected_analyses = request.form.getlist("analysis")

        print(f">>> Marker: {marker_code}")
        print(f">>> Tür listesi: {species_list}")
        print(f">>> Sekans limiti: {seq_limit}")
        print(f">>> Seçilen analizler: {selected_analyses}")

        # BOLD API'den verileri çek
        all_sequences = []
        for sp in species_list:
            try:
                print(f">>> BOLD API çağrılıyor: {sp}")
                seqs = fetch_sequences_by_organism_and_marker(
                    organism=sp,
                    marker_code=marker_code,
                    limit=seq_limit
                )
                print(f"    → {sp} için {len(seqs)} sekans bulundu.")
                all_sequences.extend(seqs)
            except Exception as e:
                print(f"    ⚠️ HATA: {sp} için veri çekilirken: {e}")

        print(f">>> Toplam çekilen sekans sayısı: {len(all_sequences)}")
        if not all_sequences:
            print("⚠️ Hiç sekans bulunamadı.")
            return render_template("align.html", error="No sequences found for the selected species and marker.")

        # FASTA dosyası oluştur
        os.makedirs(TEMP_FASTA_DIR, exist_ok=True)
        temp_fasta = os.path.join(TEMP_FASTA_DIR, "user_sequences.fasta")

        try:
            print(">>> routes.py: FASTA oluşturma başlıyor")
            final_fasta = create_fasta_from_sequences(all_sequences, temp_fasta, align=True)
            print(f">>> FASTA oluşturuldu: {final_fasta}")
        except SequenceProcessingError as e:
            print(f"❌ FASTA oluşturulamadı: {e}")
            return render_template("align.html", error=f"Failed to create FASTA: {e}")

        # Analizleri çalıştır
        results = {}

        def save_result(title, image_path):
            filename = os.path.basename(image_path)
            static_path = os.path.join(STATIC_IMG_DIR, filename)
            shutil.copyfile(image_path, static_path)
            print(f"✅ {title} başarıyla kaydedildi → {filename}")
            results[title] = url_for('static', filename=f'analysis_results/{filename}')

        if "nj" in selected_analyses:
            print("🔧 NJ analizi başlatılıyor...")
            save_result("Neighbor Joining Tree", build_nj_tree(final_fasta))

        if "mp" in selected_analyses:
            print("🔧 MP analizi başlatılıyor...")
            save_result("Maximum Parsimony Tree", build_mp_tree(final_fasta))

        if "ml" in selected_analyses:
            print("🔧 ML analizi (IQ-TREE) başlatılıyor...")
            save_result("Maximum Likelihood (IQ-TREE)", run_and_draw_bootstrap_tree(final_fasta))

        if "msn" in selected_analyses:
            print("🔧 MSN analizi başlatılıyor...")
            save_result("Minimum Spanning Network", build_msn_with_distancecalculator(final_fasta))

        if "heatmap" in selected_analyses:
            print("🔧 Heatmap analizi başlatılıyor...")
            save_result("Nucleotide Diversity Heatmap", plot_nucleotide_diversity_heatmap(final_fasta))

        print("🎯 Tüm analizler tamamlandı, sonuçlar render ediliyor.")
        return render_template("results.html", result=results)

    # GET isteği
    return render_template("align.html", result=None)
