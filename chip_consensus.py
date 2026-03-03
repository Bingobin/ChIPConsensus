#!/usr/bin/env python3
"""One-command pipeline wrapper for consensus + matrix generation."""

from __future__ import annotations

import argparse
import os

from src.consensus_lib import (
    build_consensus,
    build_support_matrix,
    load_samples_from_config,
    write_consensus_bed,
    write_matrix_tsv,
)


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
    args = parser.parse_args()

    samples = load_samples_from_config(args.config)
    consensus = build_consensus(samples)
    sample_ids, sample_vectors, support_counts = build_support_matrix(samples, consensus)

    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

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

    print(f"[OK] Samples: {len(sample_ids)}")
    print(f"[OK] Consensus regions: {len(consensus)}")
    print(f"[OK] Output dir: {out_dir}")
    print(f"[OK] Wrote: {consensus_path}")
    print(f"[OK] Wrote: {matrix_path}")


if __name__ == "__main__":
    main()
