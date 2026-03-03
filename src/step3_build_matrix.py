#!/usr/bin/env python3
"""Step 3: build support matrix (count + per-sample 0/1 flags)."""

from __future__ import annotations

import argparse
import os

from consensus_lib import (
    build_consensus,
    build_support_matrix,
    load_samples_from_config,
    write_matrix_tsv,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build consensus support matrix from sample narrowPeak files."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Sample config in CSV/TSV with header: ID and PeakFile.",
    )
    parser.add_argument(
        "--out-matrix",
        default="consensus.support.tsv",
        help="Output matrix TSV path.",
    )
    args = parser.parse_args()

    samples = load_samples_from_config(args.config)
    consensus = build_consensus(samples)
    sample_ids, sample_vectors, support_counts = build_support_matrix(samples, consensus)

    out_path = os.path.abspath(args.out_matrix)
    write_matrix_tsv(
        output_path=out_path,
        consensus=consensus,
        sample_ids=sample_ids,
        sample_vectors=sample_vectors,
        support_counts=support_counts,
    )
    print(f"[OK] Consensus regions: {len(consensus)}")
    print(f"[OK] Samples: {len(sample_ids)}")
    print(f"[OK] Support matrix written: {out_path}")


if __name__ == "__main__":
    main()
