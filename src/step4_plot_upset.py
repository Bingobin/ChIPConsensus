#!/usr/bin/env python3
"""Step 4: plot an UpSet-style figure from support matrix TSV."""

from __future__ import annotations

import argparse
import os
from collections import Counter
from typing import List, Sequence, Tuple


def load_intersection_counts(matrix_path: str) -> Tuple[List[str], List[Tuple[Tuple[int, ...], int]]]:
    with open(matrix_path, "r", encoding="utf-8") as fh:
        header_line = fh.readline().strip()
        if not header_line:
            raise ValueError("Matrix file is empty.")

        header = header_line.split("\t")
        if len(header) < 5 or header[0:4] != ["chr", "start", "end", "support_count"]:
            raise ValueError(
                "Matrix header is invalid. Expected: chr, start, end, support_count, <sample IDs...>"
            )

        sample_ids = header[4:]
        if not sample_ids:
            raise ValueError("No sample columns found in matrix.")

        combo_counter: Counter[Tuple[int, ...]] = Counter()
        for line_no, raw in enumerate(fh, start=2):
            line = raw.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4 + len(sample_ids):
                raise ValueError(
                    f"Invalid matrix line at {line_no}: expected >= {4 + len(sample_ids)} columns."
                )
            try:
                combo = tuple(int(v) for v in fields[4 : 4 + len(sample_ids)])
            except ValueError as exc:
                raise ValueError(f"Invalid 0/1 flag at line {line_no}.") from exc
            if any(v not in (0, 1) for v in combo):
                raise ValueError(f"Matrix flags must be 0/1 at line {line_no}.")
            if any(combo):
                combo_counter[combo] += 1

    ranked = sorted(combo_counter.items(), key=lambda x: (-x[1], -sum(x[0]), x[0]))
    return sample_ids, ranked


def plot_upset(
    sample_ids: Sequence[str],
    combo_counts: Sequence[Tuple[Tuple[int, ...], int]],
    output_pdf: str,
    top_n: int,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required for plotting. Install it with: pip install matplotlib"
        ) from exc

    top = list(combo_counts[:top_n]) if top_n > 0 else list(combo_counts)
    if not top:
        raise ValueError("No non-empty intersections found to plot.")

    combos = [c for c, _ in top]
    counts = [n for _, n in top]

    n_sets = len(sample_ids)
    n_cols = len(combos)
    width = max(8.0, 0.6 * n_cols + 2.0)
    height = max(5.5, 0.5 * n_sets + 2.5)

    fig, (ax_bar, ax_mat) = plt.subplots(
        2,
        1,
        figsize=(width, height),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )

    xs = list(range(n_cols))
    ax_bar.bar(xs, counts, color="#3A6EA5", width=0.8)
    ax_bar.set_ylabel("Count")
    ax_bar.set_title("Intersection Sizes (UpSet)")
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    y_positions = list(range(n_sets))
    for y in y_positions:
        ax_mat.hlines(y, -0.5, n_cols - 0.5, color="#E5E7EB", linewidth=0.8, zorder=0)

    for x, combo in enumerate(combos):
        active = [i for i, v in enumerate(combo) if v == 1]
        if not active:
            continue
        if len(active) > 1:
            ax_mat.plot([x, x], [min(active), max(active)], color="#1F2937", linewidth=1.5, zorder=2)
        ax_mat.scatter([x] * len(active), active, s=36, color="#111827", zorder=3)

    ax_mat.set_yticks(y_positions)
    ax_mat.set_yticklabels(sample_ids)
    ax_mat.set_ylabel("Samples")
    ax_mat.set_xlabel("Intersections (ranked)")
    ax_mat.set_xlim(-0.5, n_cols - 0.5)
    ax_mat.set_xticks(xs)
    ax_mat.set_xticklabels([str(i + 1) for i in xs])
    ax_mat.spines["top"].set_visible(False)
    ax_mat.spines["right"].set_visible(False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(output_pdf) or ".", exist_ok=True)
    fig.savefig(output_pdf)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot an UpSet-style intersection chart from consensus support matrix."
    )
    parser.add_argument(
        "--matrix",
        required=True,
        help="Input matrix TSV from step3_build_matrix.py.",
    )
    parser.add_argument(
        "--out-pdf",
        default="consensus.upset.pdf",
        help="Output vector PDF path.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Plot top N intersections by count (0 means all).",
    )
    args = parser.parse_args()

    sample_ids, combo_counts = load_intersection_counts(args.matrix)
    out_pdf = os.path.abspath(args.out_pdf)
    plot_upset(sample_ids, combo_counts, out_pdf, args.top_n)

    print(f"[OK] Samples: {len(sample_ids)}")
    print(f"[OK] Intersections: {len(combo_counts)}")
    print(f"[OK] UpSet PDF written: {out_pdf}")


if __name__ == "__main__":
    main()
