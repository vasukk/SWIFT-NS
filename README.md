SWIFT-NS
NS Sliding Window Identification and Fast Typing
A BLAST-based command-line tool for classification of influenza A virus NS gene sequences into Allele A or Allele B, detection of putative intermediate/chimeric sequences, and assignment to Autolin-defined NS sublineages.
---
Features
BLAST-based sliding-window allele classification
Allele A / Allele B assignment
Detection of putative intermediate/chimeric sequences
Whole-sequence and window-level summaries
Allele A/B identity and identity-difference (delta) plots
Sliding-window heatmap
Sublineage assignment using a separate BLAST database of PARNAS-selected representatives
Sublineage classification and percentage identity reported in the main summary
---
Requirements
Python 3.10+
NCBI BLAST+ (`blastn`)
Biopython
NumPy
Pandas
Matplotlib
---
Installation
Clone the repository:
```bash
git clone https://github.com/vasukk/NS-SWIFT.git
cd NS-SWIFT
```
Install Python packages:
```bash
pip install -r requirements.txt
```
Ensure NCBI BLAST+ is installed and `blastn` is available in your PATH.
The tool requires three prebuilt nucleotide BLAST databases, with the following prefixes relative to the script:
```text
blast_db/NS_Allele_A
blast_db/NS_Allele_B
NS1_lineage_blast_db/NS1_level2_lineages
```
The lineage database contains 227 PARNAS-selected representatives from 23 Autolin Level-2 sublineages. These database files must be present before running the lineage-enabled version of the tool.
---
Usage
```bash
python3 ns_swift.py \
    -i input.fasta \
    -o results
```
View all options:
```bash
python3 ns_swift.py -h
```
Note: The command above assumes the lineage-enabled script has been saved in the repository as `ns_swift.py`. If it has a different filename, use that filename instead.
---
Input
A FASTA file containing one or more influenza A virus NS nucleotide sequences.
Example:
```text
>Sequence1
ATGG...
>Sequence2
ATGG...
```
Alignment gaps in input sequences are removed before classification.
---
Output
The tool generates:
`main_summary.csv` — allele classification, sublineage prediction and identity, whole-sequence allele identities, and window counts
`window_identity_results.csv` — per-window allele comparisons, including `delta_A_minus_B`
`plots/` — per-sequence allele identity and delta plots
`window_allele_identity_heatmap.png` — summary of window-level allele assignments
Main summary columns
Column	Description
`sequence`	Query sequence identifier
`best_reference_hit`	Best whole-sequence allele reference hit
`final_classification`	Allele A, Allele B, Possible_intermediate, Ambiguous, or Unclassified
`lineage_classification`	Predicted Autolin Level-2 sublineage(s)
`lineage_identity`	Corresponding BLAST percentage identity/identities
`whole_sequence_best_A_identity`	Whole-sequence identity against Allele A database
`whole_sequence_best_B_identity`	Whole-sequence identity against Allele B database
`total_windows`	Number of sliding windows
`A_windows`	Number of Allele A windows
`B_windows`	Number of Allele B windows
`Tie_windows`	Number of tied windows
`No_hit_windows`	Number of windows without a qualifying hit
---
Workflow
```text
Input NS FASTA
      |
      v
Clean sequences and remove alignment gaps
      |
      +-----------------------------+
      |                             |
      v                             v
Generate 200-bp windows      Whole-sequence BLAST
(step size: 20 bp)           against lineage database
      |                             |
      v                             v
BLAST against Allele A       Rank sublineage hits
and Allele B databases       and report prediction
      |                             |
      v                             |
Window-level assignment             |
      |                             |
      v                             |
Sequence-level allele call          |
      |                             |
      +--------------+--------------+
                     |
                     v
         CSV summaries and plots
```
---
Allele classification
Each sequence is divided into overlapping 200-bp windows with a 20-bp step. Each window is searched independently against the Allele A and Allele B BLAST databases.
Windows are classified as:
Allele A
Allele B
Tie
No hit
The current classification settings require at least 80% identity and 100% query coverage for a window-level allele assignment.
A sequence is flagged as `Possible_intermediate` when at least two windows support Allele A and at least two support Allele B. Such a flag indicates conflicting allele signals and is not, by itself, proof of recombination.
The tool also reports per-window identity differences:
```text
delta_A_minus_B = Allele A identity - Allele B identity
```
Positive values favor Allele A; negative values favor Allele B.
---
Sublineage classification
NS sublineages were defined using Autolin on a midpoint-rooted NS phylogeny, with a minimum lineage size of 25 sequences. PARNAS was applied independently within the 23 Level-2 sublineages using a phylogenetic radius of 0.05, yielding 227 representative sequences. These representatives were used to build a nucleotide BLAST database.
Each query sequence is also searched against this representative database. The lineage classification procedure:
Retains hits meeting the minimum 80% whole-query coverage.
Identifies the best representative hit for each distinct sublineage.
Ranks candidate sublineages by BLAST bit score, with percentage identity as a secondary criterion.
Reports the selected sublineage and its percentage identity in `main_summary.csv`.
Current provisional lineage reporting rule: In the threshold-enabled version, a top hit with at least 98% identity is reported as a single sublineage. Below 98%, the two highest-ranked distinct sublineages and their corresponding identities are reported when available. This threshold is under validation and may change.
Because the hits are ranked primarily by bit score, the first reported hit does not necessarily have the highest percentage identity. A final best-versus-second-best identity-margin rule has not yet been established.
---
Future directions
Validation and refinement of sublineage assignment thresholds
Expanded analysis of intermediate sequences using conserved nucleotide and protein markers
Updated reference databases
---
Citation
Citation information will be added when available.
