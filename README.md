# Bioinformatics II - Lab Portfolio

**Student Name:** Tan Hui Shan  
**Course:** Bioinformatics II  
**Technology Stack:** Python 3.14, Streamlit, Biopython, Pandas  

---

## 📌 Lab 1: Retrieval and Basic Analysis of Protein Sequences

### Objective
Retrieve protein sequence data from NCBI Entrez / UniProt databases and perform basic physicochemical analyses using Biopython's `Bio.SeqUtils.ProtParam.ProteinAnalysis`.

### Features
- **NCBI Entrez Integration:** Retrieves FASTA records using `Bio.Entrez.efetch` and parses with `Bio.SeqIO.read`.
- **Physicochemical Properties:**
  - Sequence Length (aa)
  - Molecular Weight (in Da and kDa)
  - Isoelectric Point (pI)
  - Amino Acid Composition & Percentage
- **Interactive Visualization:**
  - Composition bar chart
  - Detailed residue frequency table
  - Collapsible FASTA & raw sequence viewer
- **Quick-Fill Examples:** One-click testing with `P04637` (Tumor protein p53), `P69905` (Hemoglobin subunit alpha), and `P01308` (Insulin).
- **Multi-Lab Sidebar Navigation:** Seamlessly switch between Lab 1 and upcoming labs in the same portfolio application.

---

## 🛠️ Setup & Execution

### 1. Install Dependencies
```bash
py -m pip install -r requirements.txt
```

### 2. Run Application
```bash
py -m streamlit run "lab1-TanHuiShan.py"
```

Then open your browser at `http://localhost:8501`.

