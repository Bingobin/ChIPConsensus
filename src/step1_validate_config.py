#!/usr/bin/env python3
"""Step 1: validate config and write normalized sample manifest."""

from __future__ import annotations

import argparse
import os

from consensus_lib import load_samples_from_config, merge_by_chrom, read_bed3_like


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate sample config (ID + narrowPeak path)."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Sample config in CSV/TSV with header: ID and PeakFile.",
    )
    parser.add_argument(
        "--out-manifest",
        default="samples.manifest.tsv",
        help="Normalized manifest output (sample_id, absolute_path, peak_count).",
    )
    args = parser.parse_args()

    samples = load_samples_from_config(args.config)
    out_path = os.path.abspath(args.out_manifest)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("sample_id\tpeak_path\tpeak_count\tmerged_peak_count\n")
        for sample in samples:
            intervals = read_bed3_like(sample.peak_path)
            peak_count = sum(len(v) for v in intervals.values())
            merged = merge_by_chrom(intervals)
            merged_count = sum(len(v) for v in merged.values())
            out.write(
                f"{sample.sample_id}\t{sample.peak_path}\t{peak_count}\t{merged_count}\n"
            )

    print(f"[OK] Config valid: {len(samples)} samples")
    print(f"[OK] Manifest written: {out_path}")


if __name__ == "__main__":
    main()
