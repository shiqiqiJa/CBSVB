# CBSVB
A tool for identifying and reconstructing bridge-like complex translocations using both second- and third-generation sequencing data

## 🧠 Background

Complex chromosomal rearrangements are central to the progression of many diseases, especially cancer. While well-known models such as *chromothripsis* have been extensively studied, less attention has been paid to patterns like *bridge-like translocations* — SVs that form via chained breakpoints connected by deletions or duplications. **CBSVB** aims to fill this gap.

---

## ✨ Features

- 🚀 **Both second- and third-generation sequencing data support**: Optimized for PacBio/Nanopore data
- 🔗 **Graph-based integration**: Connects and merges SV signatures using breakpoint graphs
- 🧩 **Bridge-like translocation detection**: Identifies and annotates BLT structures automatically

---

## 🛠️ Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/your-username/CBSVB.git
cd CBSVB
pip install -r requirements.txt
