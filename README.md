# 🧬 BARCODE-ON

**BARCODE-ON** is a user-friendly, web-based DNA barcode analysis tool built for the **BIO 496 Graduation Project** at **Boğaziçi University**. It allows users to fetch, align, and analyze DNA barcode sequences with just a few clicks.

---

## 🚀 Features

- 🔍 **Fetch sequences** from the [BOLD Systems API](https://boldsystems.org/) by species and marker
- 🔗 **Multiple alignment** using ClustalW
- 🌳 **Maximum Likelihood Trees** using IQ-TREE
- 🌐 **Minimum Spanning Networks**
- 🧪 **Nucleotide Diversity Heatmaps**
- 🎨 Dark mode support
- 📁 Result downloads for FASTA and image files

---

## 📦 Installation

```bash
git clone https://github.com/alpbay-dev/Proje-Final.git
cd Proje-Final
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

---

## ⚙️ Usage

To start the Flask server locally:

```bash
python site.py
```

Then navigate to:

```
http://127.0.0.1:5000
```

---

## 📂 Project Structure

```
Proje-Final/
│
├── app/
│   ├── routes.py               # Main Flask routes
│   ├── templates/              # HTML templates (Jinja2)
│   ├── static/                 # CSS, images, results
│   └── helpers/
│       ├── Duru/               # API & sequence utilities
│       └── Tuana/              # Analysis functions (ML, MSN, Heatmap)
│
├── site.py                    # App entry point
├── requirements.txt
└── README.md
```

---

## 👩‍💻 Team

- **Kutlu Alp Bayır** — Backend integration, sequence handling, BOLD API support  
- **Tuana Doğan** — Data analysis, tree/heatmap plotting
- **Duru Uluyol** — Sequence retrieval, FASTA creation, alignment pipeline

---

## 🎯 Motivation

We aimed to create an intuitive platform to perform commonly used bioinformatics workflows without command-line tools. This helps bridge accessibility gaps in biodiversity and genetics research.

---

## 🔗 Live Demo & GitHub

GitHub Repository: [https://github.com/alpbay-dev/Proje-Final](https://github.com/alpbay-dev/Proje-Final)

---

## 📄 License

This project is intended for academic purposes and is released without a specific license. If you wish to reuse or extend it, please credit the authors.