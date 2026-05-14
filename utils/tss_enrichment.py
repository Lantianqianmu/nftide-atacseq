#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# python tss_enrichment.py -i /data/xrz/LXYdata/out/Hep3B-1/fragments/Hep3B-1_hg38_fragments.tsv -t hg38_ArchR_TSS.tsv -o Hep3B-1_hg38 -p 32
# python tss_enrichment.py -i /data/xrz/LXYdata/out/AD38wt-adddox-48h-IFN-1/fragments/AD38wt-adddox-48h-IFN-1_hg38_fragments.tsv -t hg38_ArchR_TSS.tsv -o AD38wt-adddox-48h-IFN-1 -p 32
# python tss_enrichment.py -i /data/xrz/charseq/DS901_ATAC/PLC5-1/PLC5-1_hg38_fragments.tsv -t hg38_ArchR_TSS.tsv -o PLC5-1_hg38 -p 32


"""
bulk_atac_tss_enrichment_multiprocessing.py

Compute bulk ATAC-seq TSS enrichment from:
    1. Tn5-offset fragment TSV
       format:
           chr start end

       IMPORTANT:
           start and end are already corrected Tn5 insertion sites.
           BOTH ends will be counted independently.

    2. TSS BED file

Outputs:
    1. position-TSS enrichment txt
    2. TSS enrichment PNG

Features:
    - multiprocessing over TSS
    - strand-aware normalization
    - binary-search acceleration
    - insertion-based ATAC enrichment

"""

import argparse
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# argparse
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description="Bulk ATAC-seq TSS enrichment"
    )

    parser.add_argument(
        "-i", "--insertions",
        required=True,
        help="Tn5-offset fragment TSV: chr start end"
    )

    parser.add_argument(
        "-t", "--tss",
        required=True,
        help="TSS BED file"
    )

    parser.add_argument(
        "-o", "--outprefix",
        required=True,
        help="Output prefix"
    )

    parser.add_argument(
        "-u", "--upstream",
        type=int,
        default=2000,
        help="Upstream distance (default: 2000)"
    )

    parser.add_argument(
        "-d", "--downstream",
        type=int,
        default=2000,
        help="Downstream distance (default: 2000)"
    )

    parser.add_argument(
        "--flank",
        type=int,
        default=100,
        help="Flanking bases for normalization (default: 100)"
    )

    parser.add_argument(
        "-p", "--threads",
        type=int,
        default=4,
        help="Number of processes (default: 4)"
    )

    return parser.parse_args()


# ============================================================
# Load insertions
# ============================================================

def load_insertions(tsv_file):

    """
    Load insertion positions.

    Input TSV:
        chr start end

    start and end are both insertion sites. Note: bed files are 0-based, so minus end pos by 1.

    Returns:
        dict[chrom] -> sorted numpy array
    """

    insertions = defaultdict(list)

    with open(tsv_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split()
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # BOTH insertion sites are counted
            insertions[chrom].append(start)
            insertions[chrom].append(end - 1)

    # sort and convert to numpy
    for chrom in insertions:

        insertions[chrom] = np.array(
            sorted(insertions[chrom]),
            dtype=np.int32
        )

    return insertions


# ============================================================
# Load TSS BED
# ============================================================

def load_tss_bed(bed_file):

    """
    TSS extracted from gtf file. Note: gtf is 1-based, so minus all coordinates by 1.
    TSS structure: chr, TSS_pos, strand

    Returns:
        list of:
            (chrom, tss, strand)
    """

    tss_list = []

    with open(bed_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split()
            chrom = fields[0]
            pos = int(fields[1])
            strand = str(fields[2])
            tss = pos - 1
            tss_list.append((chrom, tss, strand))
    return tss_list


# ============================================================
# Global variables for multiprocessing
# ============================================================

GLOBAL_INSERTIONS = None
GLOBAL_UPSTREAM = None
GLOBAL_DOWNSTREAM = None
GLOBAL_WINDOW = None


def init_worker(
    insertions,
    upstream,
    downstream
):

    global GLOBAL_INSERTIONS
    global GLOBAL_UPSTREAM
    global GLOBAL_DOWNSTREAM
    global GLOBAL_WINDOW

    GLOBAL_INSERTIONS = insertions
    GLOBAL_UPSTREAM = upstream
    GLOBAL_DOWNSTREAM = downstream

    GLOBAL_WINDOW = (
        upstream +
        downstream +
        1
    )


# ============================================================
# Worker function
# ============================================================

def process_tss_chunk(tss_chunk):

    """
    Process one chunk of TSS intervals.

    Returns:
        profile
        valid_tss_count
    """

    profile = np.zeros(
        GLOBAL_WINDOW,
        dtype=np.float64
    )

    valid_tss = 0

    for chrom, tss, strand in tss_chunk:

        if chrom not in GLOBAL_INSERTIONS:
            continue

        positions = GLOBAL_INSERTIONS[chrom]

        # TSS window
        start = tss - GLOBAL_UPSTREAM
        end = tss + GLOBAL_DOWNSTREAM

        # binary search
        left_idx = np.searchsorted(
            positions,
            start,
            side="left"
        )

        right_idx = np.searchsorted(
            positions,
            end,
            side="right"
        )

        # nearby insertions
        nearby = positions[left_idx:right_idx]

        valid_tss += 1

        # accumulate insertions
        for pos in nearby:

            rel = pos - tss

            # strand normalization
            if strand == "-":
                rel = -rel

            idx = rel + GLOBAL_UPSTREAM

            if 0 <= idx < GLOBAL_WINDOW:
                profile[idx] += 1

    return profile, valid_tss


# ============================================================
# Split TSS into chunks
# ============================================================

def chunkify(lst, n):

    chunk_size = max(
        1,
        len(lst) // n
    )

    chunks = []

    for i in range(
        0,
        len(lst),
        chunk_size
    ):

        chunks.append(
            lst[i:i + chunk_size]
        )

    return chunks


# ============================================================
# Normalize profile
# ============================================================

def normalize_profile(
    profile,
    flank
):

    """
    Normalize by outer flanking regions.
    """

    left_bg = profile[:flank]

    right_bg = profile[-flank:]

    background = np.mean(
        np.concatenate([
            left_bg,
            right_bg
        ])
    )

    if background == 0:
        raise RuntimeError(
            "Background signal is zero"
        )

    enrichment = profile / background

    return enrichment


# ============================================================
# Save TXT
# ============================================================

def save_txt(
    positions,
    enrichment,
    outfile
):

    with open(outfile, "w") as out:

        out.write(
            "Position\tTSS_Enrichment\n"
        )

        for p, e in zip(
            positions,
            enrichment
        ):

            out.write(
                f"{p}\t{e:.6f}\n"
            )


# ============================================================
# Plot
# ============================================================

def plot_enrichment(
    positions,
    enrichment,
    outfile
):

    plt.figure(figsize=(8, 6))

    plt.plot(
        positions,
        enrichment,
        linewidth=1.5
    )

    plt.axvline(
        0,
        linestyle="--"
    )

    plt.xlabel(
        "Distance from TSS (bp)"
    )

    plt.ylabel(
        "Normalized TSS enrichment"
    )

    plt.title(
        "Bulk ATAC-seq TSS Enrichment"
    )

    plt.tight_layout()

    plt.savefig(
        outfile,
        dpi=300
    )

    plt.close()


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    print("Loading insertions...")

    insertions = load_insertions(
        args.insertions
    )

    total_insertions = sum(
        len(v)
        for v in insertions.values()
    )

    print(
        f"Loaded {total_insertions:,} insertions"
    )

    print("Loading TSS BED...")

    tss_list = load_tss_bed(
        args.tss
    )

    print(
        f"Loaded {len(tss_list):,} TSS"
    )

    print("Splitting TSS into chunks...")

    tss_chunks = chunkify(
        tss_list,
        args.threads
    )

    print(
        f"Running with "
        f"{args.threads} processes..."
    )

    with Pool(
        processes=args.threads,
        initializer=init_worker,
        initargs=(
            insertions,
            args.upstream,
            args.downstream
        )
    ) as pool:

        results = pool.map(
            process_tss_chunk,
            tss_chunks
        )

    print("Merging results...")

    window_size = (
        args.upstream +
        args.downstream +
        1
    )

    profile = np.zeros(
        window_size,
        dtype=np.float64
    )

    total_valid_tss = 0

    for sub_profile, sub_n in results:

        profile += sub_profile

        total_valid_tss += sub_n

    if total_valid_tss == 0:

        raise RuntimeError(
            "No valid TSS found"
        )

    # average across TSS
    profile /= total_valid_tss

    print("Normalizing profile...")

    enrichment = normalize_profile(
        profile,
        args.flank
    )

    positions = np.arange(
        -args.upstream,
        args.downstream + 1
    )

    txt_out = (
        f"{args.outprefix}"
        f"_tss_enrichment.txt"
    )

    png_out = (
        f"{args.outprefix}"
        f"_tss_enrichment.png"
    )

    print("Saving TXT...")

    save_txt(
        positions,
        enrichment,
        txt_out
    )

    print("Saving PNG...")

    plot_enrichment(
        positions,
        enrichment,
        png_out
    )

    print("Done.")

    print(f"TXT : {txt_out}")
    print(f"PNG : {png_out}")


if __name__ == "__main__":
    main()