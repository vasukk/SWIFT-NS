#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import pandas as pd
import numpy as np
import argparse
import subprocess
import tempfile
import shutil
import textwrap
import os
import re



#INPUT FILES
parser = argparse.ArgumentParser(
    description="BLAST-based sliding-window classifier for Influenza A NS alleles"
)

parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="Input FASTA file"
)

parser.add_argument(
    "-o",
    "--outdir",
    default="NS_sliding_window_results",
    help="Output directory"
)

args = parser.parse_args()

query_fasta = args.input

script_dir = os.path.dirname(os.path.abspath(__file__))

allele_a_db = os.path.join(
    script_dir,
    "blast_db",
    "NS_Allele_A"
)

allele_b_db = os.path.join(
    script_dir,
    "blast_db",
    "NS_Allele_B"
)

# SLIDING-WINDOW SETTINGS
window_size = 200
step_size = 20

# Minimum accepted BLAST percentage identity
minimum_identity = 80.0

# Require the entire 200-nt window to align
minimum_window_query_coverage = 100.0

# Require at least 80% coverage for whole-sequence BLAST
minimum_whole_query_coverage = 80.0

# Minimum identity difference required for a direct A/B call
identity_difference_threshold = 1.0

# Number of BLAST threads
num_threads = 8

# FINAL CLASSIFICATION SETTINGS


# At least 2 A windows and at least 2 B windows:
# Possible_inermediate
minimum_a_windows_for_recombinant = 2
minimum_b_windows_for_recombinant = 2


# PLOT SETTINGS
# Wrapping width for full sequence names on individual plots
plot_title_wrap_width = 100

# Wrapping width for full sequence names on heatmap rows
heatmap_label_wrap_width = 70

# Maximum heatmap dimensions in inches.
# These limits prevent very large images from exhausting memory.
maximum_heatmap_width = 40
maximum_heatmap_height = 50

# Heatmap output resolution
heatmap_dpi = 200

# OUTPUT FILES

outdir = args.outdir
plot_dir = os.path.join(outdir, "plots")

main_summary_csv = os.path.join(
    outdir,
    "main_summary.csv"
)

window_results_csv = os.path.join(
    outdir,
    "window_identity_results.csv"
)

heatmap_png = os.path.join(
    outdir,
    "window_allele_identity_heatmap.png"
)

os.makedirs(outdir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

# CHECK INPUTS

if not os.path.exists(query_fasta):
    raise FileNotFoundError(
        f"Query FASTA not found: {query_fasta}"
    )

if shutil.which("blastn") is None:
    raise RuntimeError(
        "blastn was not found. Activate the environment "
        "containing NCBI BLAST+."
    )


def blast_database_exists(prefix):

    extensions = [
        ".ndb",
        ".nhr",
        ".nin",
        ".nog",
        ".nos",
        ".not",
        ".nsq",
        ".ntf",
        ".nto"
    ]

    return any(
        os.path.exists(prefix + extension)
        for extension in extensions
    )


if not blast_database_exists(allele_a_db):
    raise FileNotFoundError(
        f"Allele A BLAST database not found: {allele_a_db}"
    )

if not blast_database_exists(allele_b_db):
    raise FileNotFoundError(
        f"Allele B BLAST database not found: {allele_b_db}"
    )



# READ AND CLEAN QUERY FASTA
query_records = []
seen_ids = set()

# Stores complete FASTA headers
description_lookup = {}


for record in SeqIO.parse(query_fasta, "fasta"):

    if record.id in seen_ids:
        print("Skipping duplicate sequence ID:", record.id)
        continue

    # Remove MAFFT gaps
    sequence = (
        str(record.seq)
        .upper()
        .replace("-", "")
    )

    # Replace unexpected ambiguity symbols with N
    sequence = re.sub(
        r"[^ACGTN]",
        "N",
        sequence
    )

    if not sequence:
        print("Skipping empty sequence:", record.id)
        continue

    cleaned_record = SeqRecord(
        Seq(sequence),
        id=record.id,
        description=record.description
    )

    query_records.append(cleaned_record)

    # record.description contains the complete FASTA header
    description_lookup[record.id] = record.description

    seen_ids.add(record.id)


if not query_records:
    raise ValueError(
        "No usable sequences were found in the query FASTA."
    )


print("Query sequences loaded:", len(query_records))

# TEMPORARY WORKING DIRECTORY
# Intermediate FASTA and BLAST files are deleted automatically

with tempfile.TemporaryDirectory(
    prefix="NS_allele_classifier_"
) as tempdir:

    cleaned_query_fasta = os.path.join(
        tempdir,
        "cleaned_queries.fasta"
    )

    window_fasta = os.path.join(
        tempdir,
        "windows.fasta"
    )

    whole_a_output = os.path.join(
        tempdir,
        "whole_sequences_vs_A.tsv"
    )

    whole_b_output = os.path.join(
        tempdir,
        "whole_sequences_vs_B.tsv"
    )

    window_a_output = os.path.join(
        tempdir,
        "windows_vs_A.tsv"
    )

    window_b_output = os.path.join(
        tempdir,
        "windows_vs_B.tsv"
    )

    SeqIO.write(
        query_records,
        cleaned_query_fasta,
        "fasta"
    )

    # BLAST FUNCTION

    blast_outfmt = (
        "6 qseqid sseqid stitle pident length qlen "
        "evalue bitscore"
    )


    def run_blast(
        query_file,
        database,
        output_file
    ):

        command = [
            "blastn",
            "-task", "blastn",
            "-query", query_file,
            "-db", database,
            "-outfmt", blast_outfmt,
            "-max_target_seqs", "20",
            "-evalue", "1e-10",
            "-num_threads", str(num_threads),
            "-out", output_file
        ]

        print("\nRunning:")
        print(" ".join(command))

        subprocess.run(
            command,
            check=True
        )

    # BLAST WHOLE SEQUENCES AGAINST BOTH DATABASES
    

    print("\n**********************************************")
    print("Blasting Whole Sequences")
    print("**********************************************")

    run_blast(
        cleaned_query_fasta,
        allele_a_db,
        whole_a_output
    )

    run_blast(
        cleaned_query_fasta,
        allele_b_db,
        whole_b_output
    )


    # READ BEST BLAST HIT

    blast_columns = [
        "query_id",
        "hit",
        "hit_title",
        "identity",
        "alignment_length",
        "query_length",
        "evalue",
        "bitscore"
    ]


    def get_best_hits(
        filename,
        prefix,
        minimum_coverage
    ):

        output_columns = [
            "query_id",
            f"{prefix}_identity",
            f"{prefix}_bitscore",
            f"{prefix}_hit",
            f"{prefix}_hit_title"
        ]

        if (
            not os.path.exists(filename)
            or os.path.getsize(filename) == 0
        ):
            return pd.DataFrame(
                columns=output_columns
            )

        blast_df = pd.read_csv(
            filename,
            sep="\t",
            names=blast_columns
        )

        numeric_columns = [
            "identity",
            "alignment_length",
            "query_length",
            "evalue",
            "bitscore"
        ]

        for column in numeric_columns:

            blast_df[column] = pd.to_numeric(
                blast_df[column],
                errors="coerce"
            )

        blast_df["query_coverage"] = (
            100.0
            * blast_df["alignment_length"]
            / blast_df["query_length"]
        )

        # Retain only hits passing both thresholds
        blast_df = blast_df[
            (
                blast_df["identity"]
                >= minimum_identity
            )
            &
            (
                blast_df["query_coverage"]
                >= minimum_coverage
            )
        ].copy()

        if blast_df.empty:
            return pd.DataFrame(
                columns=output_columns
            )

        # Select best hit:
        # 1. highest bit score
        # 2. highest identity
        # 3. longest alignment
        # 4. lowest E-value
        blast_df = blast_df.sort_values(
            by=[
                "query_id",
                "bitscore",
                "identity",
                "alignment_length",
                "evalue"
            ],
            ascending=[
                True,
                False,
                False,
                False,
                True
            ]
        )

        best_df = (
            blast_df
            .groupby(
                "query_id",
                as_index=False
            )
            .first()
        )

        best_df = best_df[[
            "query_id",
            "identity",
            "bitscore",
            "hit",
            "hit_title"
        ]]

        best_df = best_df.rename(columns={
            "identity":
                f"{prefix}_identity",

            "bitscore":
                f"{prefix}_bitscore",

            "hit":
                f"{prefix}_hit",

            "hit_title":
                f"{prefix}_hit_title"
        })

        return best_df


    
    # WHOLE-SEQUENCE BEST HITS
    
    whole_best_a = get_best_hits(
        whole_a_output,
        "whole_A",
        minimum_whole_query_coverage
    )

    whole_best_b = get_best_hits(
        whole_b_output,
        "whole_B",
        minimum_whole_query_coverage
    )


    whole_results = pd.DataFrame({
        "sequence": [
            record.id
            for record in query_records
        ],

        "full_sequence_name": [
            description_lookup[record.id]
            for record in query_records
        ]
    })


    whole_results = whole_results.merge(
        whole_best_a.rename(
            columns={
                "query_id": "sequence"
            }
        ),
        on="sequence",
        how="left"
    )


    whole_results = whole_results.merge(
        whole_best_b.rename(
            columns={
                "query_id": "sequence"
            }
        ),
        on="sequence",
        how="left"
    )


   
    # CHOOSE THE OVERALL BEST REFERENCE: ACCESSION + FULL NAME
    
    def clean_reference_accession(hit):

        if pd.isna(hit):
            return None

        hit = str(hit).strip()

        # Examples:
        # gb|GU086077.1|   -> GU086077.1
        # ref|NC_123.1|    -> NC_123.1
        parts = hit.split("|")

        if len(parts) >= 2 and parts[1]:
            return parts[1]

        return hit.rstrip("|")


    def accession_and_full_name(hit, title):

        accession = clean_reference_accession(hit)

        if accession is None and pd.isna(title):
            return "Below_threshold"

        if pd.isna(title):
            return accession

        title = str(title).strip()

        # Avoid printing the accession twice when stitle already begins
        # with forms such as gb|ACCESSION| or ACCESSION.
        removable_prefixes = [
            f"gb|{accession}|",
            f"ref|{accession}|",
            f"emb|{accession}|",
            f"dbj|{accession}|",
            accession
        ]

        for prefix in removable_prefixes:
            if prefix and title.startswith(prefix):
                title = title[len(prefix):].lstrip(" |")
                break

        if not title:
            return accession

        return f"{accession} | {title}"


    def choose_best_reference_hit(row):

        a_identity = row.get("whole_A_identity")
        b_identity = row.get("whole_B_identity")

        a_bitscore = row.get("whole_A_bitscore")
        b_bitscore = row.get("whole_B_bitscore")

        a_hit = row.get("whole_A_hit")
        b_hit = row.get("whole_B_hit")

        a_title = row.get("whole_A_hit_title")
        b_title = row.get("whole_B_hit_title")

        a_valid = pd.notna(a_identity)
        b_valid = pd.notna(b_identity)

        if not a_valid and not b_valid:
            return "Below_threshold"

        if a_valid and not b_valid:
            return accession_and_full_name(a_hit, a_title)

        if b_valid and not a_valid:
            return accession_and_full_name(b_hit, b_title)

        if a_identity > b_identity:
            return accession_and_full_name(a_hit, a_title)

        if b_identity > a_identity:
            return accession_and_full_name(b_hit, b_title)

        if pd.notna(a_bitscore) and pd.notna(b_bitscore):
            if a_bitscore > b_bitscore:
                return accession_and_full_name(a_hit, a_title)
            if b_bitscore > a_bitscore:
                return accession_and_full_name(b_hit, b_title)

        # Complete tie: report both references.
        a_reference = accession_and_full_name(a_hit, a_title)
        b_reference = accession_and_full_name(b_hit, b_title)

        if a_reference == b_reference:
            return a_reference

        return f"{a_reference} || {b_reference}"


    whole_results["best_reference_hit"] = (
        whole_results.apply(
            choose_best_reference_hit,
            axis=1
        )
    )


        # CREATE SLIDING WINDOWS
    

    print("\n**********************************************")
    print("Creating Sliding Windows")
    print("**********************************************")

    window_records = []
    window_metadata = []
    sequences_without_windows = set()


    for record in query_records:

        sequence = str(record.seq)
        sequence_length = len(sequence)

        if sequence_length < window_size:

            print(
                f"Sequence too short for windows: "
                f"{record.id} "
                f"(length={sequence_length})"
            )

            sequences_without_windows.add(
                record.id
            )

            continue


        starts = list(
            range(
                0,
                sequence_length
                - window_size
                + 1,
                step_size
            )
        )


        # Ensure the final nucleotide is covered
        final_start = (
            sequence_length
            - window_size
        )

        if final_start not in starts:
            starts.append(final_start)

        starts = sorted(set(starts))


        for window_number, start0 in enumerate(
            starts,
            start=1
        ):

            end0 = start0 + window_size

            window_sequence = sequence[
                start0:end0
            ]

            start_position = start0 + 1
            end_position = end0

            midpoint = (
                start_position
                + end_position
            ) / 2

            window_id = (
                f"{record.id}"
                f"__w{window_number}"
                f"__s{start_position}"
                f"__e{end_position}"
            )

            window_records.append(
                SeqRecord(
                    Seq(window_sequence),
                    id=window_id,
                    description=""
                )
            )

            window_metadata.append({
                "window_id": window_id,
                "sequence": record.id,
                "full_sequence_name":
                    description_lookup[record.id],
                "window": window_number,
                "start": start_position,
                "end": end_position,
                "midpoint": midpoint
            })


    if not window_records:
        raise ValueError(
            "No sliding windows were generated. "
            "Check sequence lengths or reduce window_size."
        )


    SeqIO.write(
        window_records,
        window_fasta,
        "fasta"
    )


    window_metadata_df = pd.DataFrame(
        window_metadata
    )


    print(
        "Sliding windows generated:",
        len(window_records)
    )


    
    # BLAST WINDOWS AGAINST BOTH DATABASES
    

    print("\n**********************************************")
    print("Blasting Sliding Windows")
    print("**********************************************")

    run_blast(
        window_fasta,
        allele_a_db,
        window_a_output
    )

    run_blast(
        window_fasta,
        allele_b_db,
        window_b_output
    )


    window_best_a = get_best_hits(
        window_a_output,
        "A",
        minimum_window_query_coverage
    )

    window_best_b = get_best_hits(
        window_b_output,
        "B",
        minimum_window_query_coverage
    )


    window_best_a = window_best_a.rename(
        columns={
            "query_id": "window_id"
        }
    )

    window_best_b = window_best_b.rename(
        columns={
            "query_id": "window_id"
        }
    )


    
    # COMBINE WINDOW RESULTS
    

    window_results = window_metadata_df.merge(
        window_best_a,
        on="window_id",
        how="left"
    )


    window_results = window_results.merge(
        window_best_b,
        on="window_id",
        how="left"
    )


    for column in [
        "A_identity",
        "B_identity",
        "A_bitscore",
        "B_bitscore"
    ]:

        window_results[column] = pd.to_numeric(
            window_results[column],
            errors="coerce"
        )


    
    # CLASSIFY EACH WINDOW
    

    def classify_window(row):

        a_identity = row["A_identity"]
        b_identity = row["B_identity"]

        # No qualifying hit in either database
        if (
            pd.isna(a_identity)
            and pd.isna(b_identity)
        ):
            return "No_hit"

        # Only A has a qualifying hit
        if (
            pd.notna(a_identity)
            and pd.isna(b_identity)
        ):
            return "A"

        # Only B has a qualifying hit
        if (
            pd.isna(a_identity)
            and pd.notna(b_identity)
        ):
            return "B"

        identity_difference = (
            a_identity
            - b_identity
        )

        if (
            identity_difference
            >= identity_difference_threshold
        ):
            return "A"

        if (
            identity_difference
            <= -identity_difference_threshold
        ):
            return "B"

        # Difference is below 1 percentage point.
        # Use bit score as a tie-breaker.
        a_bitscore = row["A_bitscore"]
        b_bitscore = row["B_bitscore"]

        if (
            pd.notna(a_bitscore)
            and pd.notna(b_bitscore)
        ):

            if a_bitscore > b_bitscore:
                return "A"

            if b_bitscore > a_bitscore:
                return "B"

        return "Tie"


    window_results["window_call"] = (
        window_results.apply(
            classify_window,
            axis=1
        )
    )


    window_results = window_results.sort_values(
        [
            "sequence",
            "window"
        ]
    )


    
    # SAVE WINDOW-LEVEL CSV
    

    # Keep numeric identity columns internally for classification,
    # plots, and heatmap. Only the exported CSV uses text labels.
    window_output = window_results[[
        "sequence",
        "full_sequence_name",
        "window",
        "start",
        "end",
        "midpoint",
        "A_identity",
        "B_identity",
        "window_call"
    ]].copy()

    window_output["A_identity"] = window_output["A_identity"].apply(
        lambda value: "Below_threshold"
        if pd.isna(value)
        else round(float(value), 3)
    )

    window_output["B_identity"] = window_output["B_identity"].apply(
        lambda value: "Below_threshold"
        if pd.isna(value)
        else round(float(value), 3)
    )

    window_output.to_csv(
        window_results_csv,
        index=False
    )


    
    # CREATE MAIN SUMMARY FILE
    

    summary_rows = []


    for sequence_id, group in window_results.groupby(
        "sequence"
    ):

        group = group.sort_values(
            "window"
        )

        total_windows = len(group)

        a_windows = int(
            (
                group["window_call"]
                == "A"
            ).sum()
        )

        b_windows = int(
            (
                group["window_call"]
                == "B"
            ).sum()
        )

        tie_windows = int(
            (
                group["window_call"]
                == "Tie"
            ).sum()
        )

        no_hit_windows = int(
            (
                group["window_call"]
                == "No_hit"
            ).sum()
        )


        if (
            a_windows
            >= minimum_a_windows_for_recombinant
            and b_windows
            >= minimum_b_windows_for_recombinant
        ):

            final_classification = (
                "Possible_intermediate"
            )

        elif a_windows > b_windows:

            final_classification = "Allele_A"

        elif b_windows > a_windows:

            final_classification = "Allele_B"

        elif (
            a_windows == 0
            and b_windows == 0
        ):

            final_classification = "Unclassified"

        else:

            final_classification = "Ambiguous"


        summary_rows.append({
            "sequence":
                sequence_id,

            "full_sequence_name":
                description_lookup[sequence_id],

            "final_classification":
                final_classification,

            "total_windows":
                total_windows,

            "A_windows":
                a_windows,

            "B_windows":
                b_windows,

            "Tie_windows":
                tie_windows,

            "No_hit_windows":
                no_hit_windows
        })


    # Add sequences too short to produce windows
    for sequence_id in sorted(
        sequences_without_windows
    ):

        summary_rows.append({
            "sequence":
                sequence_id,

            "full_sequence_name":
                description_lookup[sequence_id],

            "final_classification":
                "Unclassified",

            "total_windows":
                0,

            "A_windows":
                0,

            "B_windows":
                0,

            "Tie_windows":
                0,

            "No_hit_windows":
                0
        })


    summary_df = pd.DataFrame(
        summary_rows
    )


    # Add whole-sequence identities
    summary_df = summary_df.merge(
        whole_results[[
            "sequence",
            "full_sequence_name",
            "best_reference_hit",
            "whole_A_identity",
            "whole_B_identity"
        ]],
        on=[
            "sequence",
            "full_sequence_name"
        ],
        how="outer"
    )


    summary_df = summary_df.rename(columns={
        "whole_A_identity":
            "whole_sequence_best_A_identity",

        "whole_B_identity":
            "whole_sequence_best_B_identity"
    })


    summary_df[
        "final_classification"
    ] = summary_df[
        "final_classification"
    ].fillna(
        "Unclassified"
    )


    for column in [
        "total_windows",
        "A_windows",
        "B_windows",
        "Tie_windows",
        "No_hit_windows"
    ]:

        summary_df[column] = (
            summary_df[column]
            .fillna(0)
            .astype(int)
        )


    summary_df = summary_df[[
        "sequence",
        "best_reference_hit",
        "final_classification",
        "whole_sequence_best_A_identity",
        "whole_sequence_best_B_identity",
        "total_windows",
        "A_windows",
        "B_windows",
        "Tie_windows",
        "No_hit_windows"
    ]]


    summary_df = summary_df.sort_values(
        [
            "final_classification",
            "sequence"
        ]
    )

    # Export readable text instead of blank cells when no whole-sequence
    # hit passed the identity and coverage thresholds.
    summary_output = summary_df.copy()

    summary_output["best_reference_hit"] = (
        summary_output["best_reference_hit"]
        .fillna("Below_threshold")
    )


    summary_output["whole_sequence_best_A_identity"] = (
        summary_output["whole_sequence_best_A_identity"].apply(
            lambda value: "Below_threshold"
            if pd.isna(value)
            else round(float(value), 3)
        )
    )

    summary_output["whole_sequence_best_B_identity"] = (
        summary_output["whole_sequence_best_B_identity"].apply(
            lambda value: "Below_threshold"
            if pd.isna(value)
            else round(float(value), 3)
        )
    )

    summary_output.to_csv(
        main_summary_csv,
        index=False
    )


   
    # CREATE ONE CURVE PLOT FOR EVERY SEQUENCE
    

    classification_lookup = (
        summary_df
        .set_index("sequence")
        .to_dict("index")
    )


    for sequence_id, group in window_results.groupby(
        "sequence"
    ):

        group = group.sort_values(
            "start"
        )

        information = classification_lookup[
            sequence_id
        ]

        full_sequence_name = description_lookup[sequence_id]

        classification = information[
            "final_classification"
        ]

        a_windows = information[
            "A_windows"
        ]

        b_windows = information[
            "B_windows"
        ]

        tie_windows = information[
            "Tie_windows"
        ]

        no_hit_windows = information[
            "No_hit_windows"
        ]


        wrapped_full_name = textwrap.fill(
            str(full_sequence_name),
            width=plot_title_wrap_width
        )


        plt.figure(
            figsize=(12, 7)
        )


        plt.plot(
            group["midpoint"],
            group["A_identity"],
            marker="o",
            linewidth=2,
            markersize=4,
            label="Best Allele A identity"
        )


        plt.plot(
            group["midpoint"],
            group["B_identity"],
            marker="o",
            linewidth=2,
            markersize=4,
            label="Best Allele B identity"
        )


        plt.axhline(
            y=minimum_identity,
            linestyle="--",
            linewidth=1,
            label=(
                f"{minimum_identity:.0f}% "
                "identity threshold"
            )
        )


        plt.xlabel(
            "Window midpoint position"
        )

        plt.ylabel(
            "Best-hit identity (%)"
        )


        valid_values = pd.concat([
            group["A_identity"],
            group["B_identity"]
        ]).dropna()


        if valid_values.empty:

            plt.ylim(
                max(
                    0,
                    minimum_identity - 10
                ),
                101
            )

        else:

            lower_limit = max(
                0,
                min(
                    minimum_identity - 5,
                    np.floor(
                        valid_values.min() / 5
                    ) * 5 - 5
                )
            )

            # Always leave a little headroom above 100%
            plt.ylim(
                lower_limit,
                101
            )


        plt.grid(
            alpha=0.3
        )


        # Full sequence name is included here
        plt.title(
            f"{wrapped_full_name}\n\n"
            f"Classification: {classification} | "
            f"A: {a_windows} | "
            f"B: {b_windows} | "
            f"Tie: {tie_windows} | "
            f"No hit: {no_hit_windows}",
            fontsize=10,
            pad=18
        )


        plt.legend()
        plt.tight_layout()


        safe_name = re.sub(
            r"[^A-Za-z0-9._-]+",
            "_",
            sequence_id
        )


        plot_file = os.path.join(
            plot_dir,
            f"{safe_name}_A_vs_B.png"
        )


        plt.savefig(
            plot_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


        # CREATE HEATMAP
    

    print("\n**********************************************")
    print("Creating Window Heatmap")
    print("**********************************************")


    sequence_order = (
        summary_df[
            summary_df["total_windows"] > 0
        ]["sequence"]
        .tolist()
    )


    if sequence_order:

        maximum_windows = int(
            window_results["window"].max()
        )

        number_of_sequences = len(
            sequence_order
        )


        print(
            "Heatmap sequences:",
            number_of_sequences
        )

        print(
            "Heatmap maximum windows:",
            maximum_windows
        )


        # Use uint8 rather than float64.
        # This greatly reduces memory usage.
        heatmap_rgb = np.full(
            (
                number_of_sequences,
                maximum_windows,
                3
            ),
            255,
            dtype=np.uint8
        )


        heatmap_labels = []


        for row_index, sequence_id in enumerate(
            sequence_order
        ):

            full_name = description_lookup[
                sequence_id
            ]

            wrapped_label = textwrap.fill(
                str(full_name),
                width=heatmap_label_wrap_width
            )

            heatmap_labels.append(
                wrapped_label
            )


            sequence_windows = (
                window_results[
                    window_results["sequence"]
                    == sequence_id
                ]
                .sort_values("window")
            )


            for _, row in sequence_windows.iterrows():

                column_index = (
                    int(row["window"]) - 1
                )

                window_call = row[
                    "window_call"
                ]


                if window_call == "A":

                    identity = row[
                        "A_identity"
                    ]

                    intensity = (
                        identity
                        - minimum_identity
                    ) / (
                        100.0
                        - minimum_identity
                    )

                    intensity = float(
                        np.clip(
                            intensity,
                            0.0,
                            1.0
                        )
                    )

                    # Pale blue to dark blue
                    red_value = int(
                        255
                        - 215 * intensity
                    )

                    green_value = int(
                        255
                        - 165 * intensity
                    )

                    heatmap_rgb[
                        row_index,
                        column_index
                    ] = [
                        red_value,
                        green_value,
                        255
                    ]


                elif window_call == "B":

                    identity = row[
                        "B_identity"
                    ]

                    intensity = (
                        identity
                        - minimum_identity
                    ) / (
                        100.0
                        - minimum_identity
                    )

                    intensity = float(
                        np.clip(
                            intensity,
                            0.0,
                            1.0
                        )
                    )

                    # Pale red to dark red
                    green_value = int(
                        255
                        - 190 * intensity
                    )

                    blue_value = int(
                        255
                        - 190 * intensity
                    )

                    heatmap_rgb[
                        row_index,
                        column_index
                    ] = [
                        255,
                        green_value,
                        blue_value
                    ]


                elif window_call == "Tie":

                    heatmap_rgb[
                        row_index,
                        column_index
                    ] = [
                        155,
                        155,
                        155
                    ]


                elif window_call == "No_hit":

                    heatmap_rgb[
                        row_index,
                        column_index
                    ] = [
                        255,
                        255,
                        255
                    ]


        # Determine figure size, but cap it to prevent memory failure
        figure_width = min(
            maximum_heatmap_width,
            max(
                14,
                maximum_windows * 0.25
            )
        )

        figure_height = min(
            maximum_heatmap_height,
            max(
                7,
                number_of_sequences * 0.38
            )
        )


        fig, ax = plt.subplots(
            figsize=(
                figure_width,
                figure_height
            )
        )


        ax.imshow(
            heatmap_rgb,
            aspect="auto",
            interpolation="nearest"
        )


        ax.set_xlabel(
            "Sliding-window number",
            fontsize=12
        )

        ax.set_ylabel(
            "Full sequence name",
            fontsize=12
        )


        # Show a reasonable number of x-axis labels
        if maximum_windows <= 40:

            tick_positions = np.arange(
                maximum_windows
            )

        else:

            tick_step = max(
                1,
                int(
                    np.ceil(
                        maximum_windows / 25
                    )
                )
            )

            tick_positions = np.arange(
                0,
                maximum_windows,
                tick_step
            )


        ax.set_xticks(
            tick_positions
        )

        ax.set_xticklabels(
            tick_positions + 1,
            rotation=90,
            fontsize=8
        )


        ax.set_yticks(
            np.arange(
                number_of_sequences
            )
        )

        # Full sequence names are included on the left side
        ax.set_yticklabels(
            heatmap_labels,
            fontsize=6
        )


        # Add the final sequence classification on the right side
        classification_lookup_heatmap = (
            summary_df
            .set_index("sequence")[
                "final_classification"
            ]
            .to_dict()
        )

        heatmap_classification_labels = [
            classification_lookup_heatmap.get(
                sequence_id,
                "Unclassified"
            )
            for sequence_id in sequence_order
        ]

        classification_axis = ax.twinx()

        # Keep the right-side labels aligned with heatmap rows
        classification_axis.set_ylim(
            ax.get_ylim()
        )

        classification_axis.set_yticks(
            np.arange(
                number_of_sequences
            )
        )

        classification_axis.set_yticklabels(
            heatmap_classification_labels,
            fontsize=7
        )

        classification_axis.set_ylabel(
            "Final classification",
            fontsize=12,
            labelpad=12
        )

        classification_axis.tick_params(
            axis="y",
            length=0,
            pad=6
        )


        ax.set_title(
            "Sliding-window identity heatmap\n",
            fontsize=14,
            pad=15
        )


        legend_elements = [
            Patch(
                facecolor=(
                    40 / 255,
                    90 / 255,
                    1.0
                ),
                label="Allele A window"
            ),

            Patch(
                facecolor=(
                    1.0,
                    65 / 255,
                    65 / 255
                ),
                label="Allele B window"
            ),

            Patch(
                facecolor=(
                    155 / 255,
                    155 / 255,
                    155 / 255
                ),
                label="Tie"
            ),

            Patch(
                facecolor="white",
                edgecolor="black",
                label="No hit / no window"
            )
        ]


        ax.legend(
            handles=legend_elements,
            bbox_to_anchor=(1.35, 1.0),
            loc="upper left",
            fontsize=9,
            frameon=True
        )


        # Leave room for full names without creating a huge figure
        # Leave room for names on the left, classifications on the right,
        # and the legend beyond the classification labels.
        fig.subplots_adjust(
            left=0.40,
            right=0.62,
            bottom=0.14,
            top=0.90
        )


        plt.savefig(
            heatmap_png,
            dpi=heatmap_dpi,
            bbox_inches="tight"
        )

        plt.close(fig)


        print(
            "Heatmap saved:",
            heatmap_png
        )


    else:

        print(
            "No sequences had sliding windows. "
            "Heatmap was not created."
        )



# FINAL SUMMARY

print("\n**********************************************")
print("Final Classification Counts")
print("**********************************************")

print(
    summary_df[
        "final_classification"
    ].value_counts(
        dropna=False
    )
)


print("\nIdentity threshold:")
print(
    f"Only hits with identity >= "
    f"{minimum_identity:.1f}% are accepted."
)


print("\nClassification rule:")
print(
    "Possible_intermediate = at least "
    f"{minimum_a_windows_for_recombinant} "
    "Allele A windows and at least "
    f"{minimum_b_windows_for_recombinant} "
    "Allele B windows."
)


print("\n outputs:")

print(
    "Main summary:",
    main_summary_csv
)

print(
    "Window results:",
    window_results_csv
)

print(
    "Heatmap:",
    heatmap_png
)

print(
    "Individual sequence plots:",
    plot_dir
)

print("\nDone.")
