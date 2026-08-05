# NS-SWIFT

**NS Sliding Window Identification and Fast Typing**

A BLAST-based sliding-window command-line tool for rapid classification of Influenza A virus NS gene sequences into Allele A and Allele B and identification of putative intermediate/chimeric sequences.

---

## Features

- BLAST-based sliding-window classification
- Allele A / Allele B assignment
- Detection of putative intermediate/chimeric sequences
- Whole-sequence and window-level summaries
- Identity plots for each sequence
- Sliding-window heatmap

---

## Requirements

- Python 3.10+
- NCBI BLAST+
- Biopython
- NumPy
- Pandas
- Matplotlib

---

## Installation

Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/SWAN.git
cd SWAN
```

Install Python packages

```bash
pip install -r requirements.txt
```

Ensure NCBI BLAST+ is installed and `blastn` is available in your PATH.

---

## Usage

```bash
python3 ns_swift.py \
    -i input.fasta \
    -o results
```

View all options

```bash
python3 ns_swift.py -h
```

---

## Input

A FASTA file containing one or more Influenza A virus NS nucleotide sequences.

Example:

```text
>Sequence1
ATGG...

>Sequence2
ATGG...
```

---

## Output

The tool generates:

- `main_summary.csv`
- `window_identity_results.csv`
- `plots/`
- `window_allele_identity_heatmap.png`

---

## Workflow

```
Input FASTA
      │
      ▼
Remove alignment gaps
      │
      ▼
Generate sliding windows
      │
      ▼
BLAST against Allele A database
BLAST against Allele B database
      │
      ▼
Window classification
      │
      ▼
Sequence classification
      │
      ▼
Plots + CSV summaries
```

---

## Classification

Each sequence is divided into overlapping windows.

Each window is compared against curated Allele A and Allele B BLAST databases.

Windows are classified as:

- Allele A
- Allele B
- Tie
- No hit

Sequences containing both Allele A and Allele B windows are reported as **Possible_intermediate**.

---

## Future directions

- Expanded analysis of intermediate sequences using conserved nucleotide and protein markers
- Updated reference databases

---

## Citation
