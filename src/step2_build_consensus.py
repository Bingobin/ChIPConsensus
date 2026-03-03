#!/usr/bin/env python3
"""Step 2: build merged consensus regions."""

from __future__ import annotations

import argparse
import os

from consensus_lib import build_consensus, load_samples_from_config, write_consensus_bed


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build consensus regions from multiple narrowPeak files."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Sample config in CSV/TSV with header: ID and PeakFile.",
    )
    parser.add_argument(
        "--out-consensus",
        default="consensus.bed",
        help="Output BED path for merged consensus regions.",
    )
    args = parser.parse_args()

    samples = load_samples_from_config(args.config)
    consensus = build_consensus(samples)
    out_path = os.path.abspath(args.out_consensus)
    write_consensus_bed(consensus, out_path)
    print(f"[OK] Consensus regions: {len(consensus)}")
    print(f"[OK] Consensus BED written: {out_path}")


if __name__ == "__main__":
    main()
