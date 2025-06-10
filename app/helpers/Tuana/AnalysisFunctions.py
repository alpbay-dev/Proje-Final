from Bio import AlignIO
import os
import subprocess
from Bio import Phylo
import matplotlib.pyplot as plt

file_path = "aligned_fasta_final.fasta"

def read_alignment(file_path):
    """
    Fasta fileını okuyor
    """
    alignment = AlignIO.read(file_path, "fasta")
    return alignment


def run_and_draw_bootstrap_tree(input_path: str, iqtree_path: str = "iqtree\iqtree3.exe", bootstrap: int = 1000) -> str:
    """
    IQ-TREE ile bootstrap destekli ML ağacı oluşturur ve eksensiz PNG olarak çizer.

    Args:
        input_path (str): FASTA hizalama dosyasının tam yolu (.fas)
        iqtree_path (str): iqtree3.exe dosyasının tam yolu
        bootstrap (int): Bootstrap tekrar sayısı (default: 1000)

    Returns:
        str: Kaydedilen PNG dosyasının yolu
    """

    input_path = os.path.abspath(input_path)
    work_dir = os.path.dirname(input_path)
    input_filename = os.path.basename(input_path)  # 👈 sadece dosya adı
    contree_path = os.path.join(work_dir, input_filename + ".contree")  # .contree dosyası
    output_img = os.path.join(work_dir, input_filename + "_bootstrap_tree.png")  # çıktı görseli

    # IQ-TREE komutu
    cmd = [
        iqtree_path,
        "-s", input_path,
        "-m", "GTR+G",
        "-nt", "AUTO",
        "-bb", str(bootstrap),
        "-redo"  # 🔥 önceki çıktıları silerek yeniden çalıştırır
    ]

    print("▶ IQ-TREE komutu:", " ".join(cmd))
    print("▶ Çalışma dizini:", work_dir)

    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)

    print("▶ STDOUT:\n", result.stdout)
    print("▶ STDERR:\n", result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"IQ-TREE hata verdi! Kod: {result.returncode}\n{result.stderr}")

    if not os.path.exists(contree_path):
        raise FileNotFoundError(f".contree dosyası bulunamadı: {contree_path}")

    tree = Phylo.read(contree_path, "newick")

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.set_axis_off()

    plt.savefig(output_img, dpi=300, bbox_inches="tight")
    plt.close()

    print("✅ Ağaç başarıyla çizildi:", output_img)
    return output_img


#run_and_draw_bootstrap_tree(file_path)

def build_msn_with_distancecalculator(fasta_path: str, output_img: str = "msn_distance_network.png"):
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    import networkx as nx
    import matplotlib.pyplot as plt
    from collections import defaultdict

    alignment = AlignIO.read(fasta_path, "fasta")
    seq_to_ids = defaultdict(list)
    for record in alignment:
        seq = str(record.seq)
        seq_to_ids[seq].append(record.id)

    haplotypes = list(seq_to_ids.keys())
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    aln_len = alignment.get_alignment_length()

    G = nx.Graph()

    # Düğümleri orijinal ID ile ekle
    seq_to_rep_id = {}
    for seq, ids in seq_to_ids.items():
        rep_id = ids[0]  # ilk ID (örneğin Falpe1)
        seq_to_rep_id[seq] = rep_id
        label = f"{rep_id}\n({len(ids)}x)" if len(ids) > 1 else rep_id
        G.add_node(rep_id, label=label, size=len(ids))

    # Kenarları ekle
    for i in range(len(alignment)):
        for j in range(i + 1, len(alignment)):
            seq1 = str(alignment[i].seq)
            seq2 = str(alignment[j].seq)
            rep1 = seq_to_rep_id[seq1]
            rep2 = seq_to_rep_id[seq2]
            if rep1 == rep2:
                continue

            dist_prop = distance_matrix[alignment[i].id, alignment[j].id]
            dist_nt = round(dist_prop * aln_len)
            G.add_edge(rep1, rep2, weight=dist_nt)

    MST = nx.minimum_spanning_tree(G)
    pos = nx.kamada_kawai_layout(MST, weight="weight")
    sizes = [MST.nodes[n]["size"] * 500 for n in MST.nodes]
    labels = {n: MST.nodes[n]["label"] for n in MST.nodes}
    edge_labels = {(u, v): MST[u][v]["weight"] for u, v in MST.edges}

    plt.figure(figsize=(10, 10))
    nx.draw(MST, pos, with_labels=False, node_size=sizes, node_color="lightgreen", edge_color="gray")
    nx.draw_networkx_labels(MST, pos, labels=labels, font_size=9)
    nx.draw_networkx_edge_labels(MST, pos, edge_labels=edge_labels, font_size=8)
    plt.title("Minimum Spanning Network", fontsize=14)
    plt.axis("off")
    plt.savefig(output_img, bbox_inches="tight", dpi=300)
    plt.close()

    print(f"MSN çizildi (orijinal ID etiketli): {output_img}")
    return output_img


#build_msn_with_distancecalculator(file_path)

def plot_nucleotide_diversity_heatmap(fasta_path: str, output_img: str = None):
    """
    FASTA hizalamasındaki örnekler arası nükleotid çeşitliliğine göre heatmap çizer ve PNG'ye kaydeder.
    output_img verilmezse otomatik oluşturur.

    Args:
        fasta_path (str): MSA dosyasının yolu (FASTA formatında)
        output_img (str): PNG olarak kaydedilecek dosya (opsiyonel)

    Returns:
        str: Kaydedilen PNG dosyasının yolu
    """
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import os

    print(f" Hizalama dosyası okunuyor: {fasta_path}")
    aln = AlignIO.read(fasta_path, "fasta")

    print(" Distance matrix hesaplanıyor (π)...")
    calculator = DistanceCalculator("identity")
    dist_matrix = calculator.get_distance(aln)

    print(" DataFrame oluşturuluyor...")
    df = pd.DataFrame(dist_matrix.matrix,
                      index=dist_matrix.names,
                      columns=dist_matrix.names)

    print(" Heatmap çiziliyor...")
    plt.figure(figsize=(10, 8))
    sns.heatmap(df, cmap="viridis", annot=False, square=True,
                cbar_kws={"label": "Nucleotide diversity (π)"})
    plt.title("Pairwise Nucleotide Diversity Heatmap")
    plt.tight_layout()

    #  Eğer output verilmemişse otomatik oluştur
    if output_img is None:
        base = os.path.splitext(fasta_path)[0]
        output_img = base + "_pi_heatmap.png"

    print(f" PNG olarak kaydediliyor: {output_img}")
    plt.savefig(output_img, dpi=300, bbox_inches="tight")
    plt.close()
    print(" Başarıyla kaydedildi.")

    return output_img


