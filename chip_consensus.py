#!/usr/bin/env python3
"""One-command pipeline wrapper for consensus + matrix generation."""

from __future__ import annotations

import argparse
import os

from src.consensus_lib import (
    build_consensus,
    build_support_matrix,
    load_samples_from_config,
    merge_by_chrom,
    read_bed3_like,
    write_consensus_bed,
    write_matrix_tsv,
)
from src.step4_plot_upset import load_intersection_counts, plot_upset


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run full consensus peak pipeline in pure Python."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Sample config in CSV/TSV with header: ID and PeakFile.",
    )
    parser.add_argument(
        "--prefix",
        default="chip_consensus",
        help="Output file prefix (written under --out-dir).",
    )
    parser.add_argument(
        "--out-dir",
        default=".",
        help="Root output directory.",
    )
    parser.add_argument(
        "--upset-top-n",
        type=int,
        default=20,
        help="Top N intersections to plot in UpSet PDF (0 means all).",
    )
    parser.add_argument(
        "--no-upset",
        action="store_true",
        help="Skip UpSet plotting step.",
    )
    args = parser.parse_args()

    print("[INFO] Step 1/4: load and validate config...", flush=True)
    samples = load_samples_from_config(args.config)

    print("[INFO] Step 2/4: summarize samples and build manifest...", flush=True)
    manifest_rows = []
    for sample in samples:
        intervals = read_bed3_like(sample.peak_path)
        peak_count = sum(len(v) for v in intervals.values())
        merged = merge_by_chrom(intervals)
        merged_count = sum(len(v) for v in merged.values())
        manifest_rows.append((sample.sample_id, sample.peak_path, peak_count, merged_count))

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(out_dir, f"{args.prefix}.samples.manifest.tsv")
    with open(manifest_path, "w", encoding="utf-8") as out:
        out.write("sample_id\tpeak_path\tpeak_count\tmerged_peak_count\n")
        for sample_id, peak_path, peak_count, merged_count in manifest_rows:
            out.write(f"{sample_id}\t{peak_path}\t{peak_count}\t{merged_count}\n")

    print(
        f"[INFO] Manifest ready: {len(manifest_rows)} samples -> {manifest_path}",
        flush=True,
    )
    for sample_id, _peak_path, peak_count, merged_count in manifest_rows:
        print(
            f"[INFO]   {sample_id}: peak_count={peak_count}, merged_peak_count={merged_count}",
            flush=True,
        )

    print("[INFO] Step 3/4: build consensus and support matrix...", flush=True)
    consensus = build_consensus(samples)
    sample_ids, sample_vectors, support_counts = build_support_matrix(samples, consensus)

    consensus_path = os.path.join(out_dir, f"{args.prefix}.consensus.bed")
    matrix_path = os.path.join(out_dir, f"{args.prefix}.consensus.support.tsv")
    write_consensus_bed(consensus, consensus_path)
    write_matrix_tsv(
        output_path=matrix_path,
        consensus=consensus,
        sample_ids=sample_ids,
        sample_vectors=sample_vectors,
        support_counts=support_counts,
    )

    upset_path = os.path.join(out_dir, f"{args.prefix}.consensus.upset.pdf")
    if args.no_upset:
        print("[INFO] Step 4/4: skipped UpSet plotting (--no-upset).", flush=True)
    else:
        print("[INFO] Step 4/4: plot UpSet PDF...", flush=True)
        try:
            upset_sample_ids, combo_counts = load_intersection_counts(matrix_path)
            plot_upset(upset_sample_ids, combo_counts, upset_path, args.upset_top_n)
            print(f"[OK] UpSet PDF written: {upset_path}", flush=True)
        except Exception as exc:
            print(f"[WARN] UpSet plotting failed: {exc}", flush=True)

    print(f"[OK] Samples: {len(sample_ids)}")
    print(f"[OK] Consensus regions: {len(consensus)}")
    print(f"[OK] Output dir: {out_dir}")
    print(f"[OK] Wrote: {manifest_path}")
    print(f"[OK] Wrote: {consensus_path}")
    print(f"[OK] Wrote: {matrix_path}")


if __name__ == "__main__":
    main()
