from Bio import AlignIO

file_path = "HW1_TuanaDogan_alignment_file.fas"

def read_alignment(file_path):
    """
    Fasta fileını okuyor
    """
    alignment = AlignIO.read(file_path, "fasta")
    return alignment

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

def build_nj_tree(input_path: str) -> str:
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    import matplotlib.pyplot as plt
    import os

    # 1. Alignment oku
    aln = AlignIO.read(input_path, "fasta")

    # 2. Mesafe matrisi ve NJ ağacı
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # 3. İç düğümlerde isim olmasın
    for clade in tree.get_nonterminals():
        clade.name = ""

    # 4. Çizim için eksensiz figür oluştur
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.set_axis_off()  # 🔥 X/Y eksenlerini kaldır

    # 5. Kaydet
    output_path = input_path + "_nj_tree.png"
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()

    print(f" NJ ağacı çizildi (eksensiz): {output_path}")
    return output_path


#build_nj_tree(file_path)


def build_mp_tree(input_path: str) -> str:
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
    import matplotlib.pyplot as plt
    import os

    # 1. Alignment oku
    aln = AlignIO.read(input_path, "fasta")

    # 2. Başlangıç ağacı için mesafe matrisi ve UPGMA ağacı oluştur
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(aln)
    upgma_tree = DistanceTreeConstructor().upgma(dm)

    # 3. Parsimony ağacı oluştur
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    mp_tree = ParsimonyTreeConstructor(searcher, upgma_tree).build_tree(aln)

    # 4. İç düğüm isimlerini kaldır
    for clade in mp_tree.get_nonterminals():
        clade.name = ""

    # 5. Eksensiz çizim ve kayıt
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(mp_tree, do_show=False, axes=ax)
    ax.set_axis_off()  # 👈 Ekseni kaldır

    # 6. Kaydet
    output_path = input_path + "_mp_tree.png"
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()

    print(f" MP ağacı çizildi (eksensiz): {output_path}")
    return output_path

#build_mp_tree(file_path)

import os
import subprocess
from Bio import Phylo
import matplotlib.pyplot as plt

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
    import os
    import subprocess
    from Bio import Phylo
    import matplotlib.pyplot as plt

    # 1. Yol hazırlıkları
    input_path = os.path.abspath(input_path)
    work_dir = os.path.dirname(input_path)
    contree_path = input_path + ".contree"
    output_img = input_path + "_bootstrap_tree.png"

    # 2. IQ-TREE komutunu çalıştır
    cmd = [
        iqtree_path,
        "-s", input_path,
        "-m", "GTR+G",
        "-nt", "AUTO",
        "-bb", str(bootstrap)
    ]

    print("▶ IQ-TREE çalıştırılıyor:", " ".join(cmd))
    subprocess.run(cmd, cwd=work_dir, check=True)

    # 3. .contree dosyası oluştu mu?
    if not os.path.exists(contree_path):
        raise FileNotFoundError(f" .contree dosyası bulunamadı: {contree_path}")

    # 4. Ağaç dosyasını oku
    tree = Phylo.read(contree_path, "newick")

    # 5. Görsel olarak çiz
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)

    #  X/Y eksenlerini kaldır
    ax.set_axis_off()

    # 6. Kaydet ve kapat
    plt.savefig(output_img, dpi=300, bbox_inches="tight")
    plt.close()

    print(" Bootstrap destekli ağaç çizildi (eksensiz):", output_img)
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
        pd.DataFrame: Nükleotid çeşitliliği matrisi
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

    return df







#plot_nucleotide_diversity_heatmap(file_path)


