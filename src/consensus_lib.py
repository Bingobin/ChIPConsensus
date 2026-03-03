#!/usr/bin/env python3
"""Shared utilities for consensus peak construction."""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


Interval = Tuple[int, int]
ConsensusInterval = Tuple[str, int, int]


@dataclass(frozen=True)
class Sample:
    sample_id: str
    peak_path: str


def chrom_sort_key(chrom: str) -> Tuple[int, int, str]:
    """Sort chromosomes in a biologically reasonable order."""
    c = chrom
    if c.startswith("chr"):
        c = c[3:]

    if c.isdigit():
        return (0, int(c), chrom)

    if c == "X":
        return (1, 23, chrom)
    if c == "Y":
        return (1, 24, chrom)
    if c in {"M", "MT"}:
        return (1, 25, chrom)

    m = re.match(r"^(\d+)_", c)
    if m:
        return (2, int(m.group(1)), chrom)

    return (3, 0, chrom)


def load_samples_from_config(config_path: str) -> List[Sample]:
    samples: List[Sample] = []
    seen_ids = set()
    delimiter: str | None = None
    header_checked = False

    with open(config_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "," in line:
                delimiter = ","
            elif "\t" in line:
                delimiter = "\t"
            else:
                raise ValueError(
                    "Invalid config delimiter. Use CSV (',') or TSV ('\\t') with header: ID,PeakFile."
                )
            break

    if delimiter is None:
        raise ValueError("No valid header found. First non-comment line must be: ID,PeakFile.")

    with open(config_path, "r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            if delimiter == ",":
                fields = next(csv.reader([line], skipinitialspace=True))
                fields = [f.strip() for f in fields]
            else:
                fields = [f.strip() for f in line.split("\t")]

            if len(fields) != 2:
                raise ValueError(
                    f"Invalid config at line {line_no}: expected 2 columns, got {len(fields)}"
                )

            if not header_checked:
                if fields[0].lower() != "id" or fields[1].lower() != "peakfile":
                    raise ValueError(
                        "Config header must be exactly: ID,PeakFile (CSV) or ID<TAB>PeakFile (TSV)."
                    )
                header_checked = True
                continue

            sample_id = fields[0]
            peak_path = fields[1]
            abs_peak_path = os.path.abspath(peak_path)
            if sample_id in seen_ids:
                raise ValueError(f"Duplicate sample ID '{sample_id}' at line {line_no}")
            if not os.path.exists(abs_peak_path):
                raise FileNotFoundError(
                    f"Peak file not found for sample '{sample_id}': {abs_peak_path}"
                )
            seen_ids.add(sample_id)
            samples.append(Sample(sample_id=sample_id, peak_path=abs_peak_path))

    if not header_checked:
        raise ValueError("No valid header found. First non-comment line must be: ID,PeakFile.")
    if not samples:
        raise ValueError("No samples found in config.")
    return samples


def read_bed3_like(path: str) -> Dict[str, List[Interval]]:
    """Read first three columns from BED/narrowPeak-like files."""
    chrom_to_intervals: Dict[str, List[Interval]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"Invalid interval line in {path}:{line_no}, expected >=3 columns."
                )
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid start/end in {path}:{line_no}: {fields[1]} {fields[2]}"
                ) from exc
            if start < 0 or end < 0 or end <= start:
                continue
            chrom_to_intervals.setdefault(chrom, []).append((start, end))
    return chrom_to_intervals


def merge_intervals(intervals: Sequence[Interval]) -> List[Interval]:
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged: List[Interval] = [sorted_intervals[0]]
    for start, end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def merge_by_chrom(chrom_to_intervals: Dict[str, List[Interval]]) -> Dict[str, List[Interval]]:
    return {chrom: merge_intervals(v) for chrom, v in chrom_to_intervals.items()}


def build_consensus(samples: Sequence[Sample]) -> List[ConsensusInterval]:
    union_by_chrom: Dict[str, List[Interval]] = {}
    for sample in samples:
        peaks = read_bed3_like(sample.peak_path)
        for chrom, intervals in peaks.items():
            union_by_chrom.setdefault(chrom, []).extend(intervals)

    consensus: List[ConsensusInterval] = []
    for chrom in sorted(union_by_chrom.keys(), key=chrom_sort_key):
        merged = merge_intervals(union_by_chrom[chrom])
        consensus.extend((chrom, s, e) for s, e in merged)
    return consensus


def intervals_support_vector(
    consensus_intervals: Sequence[Interval], sample_intervals: Sequence[Interval]
) -> List[int]:
    """Return 0/1 support vector for one chromosome using two-pointer sweep."""
    vector = [0] * len(consensus_intervals)
    if not consensus_intervals or not sample_intervals:
        return vector

    j = 0
    n = len(sample_intervals)
    for i, (c_start, c_end) in enumerate(consensus_intervals):
        while j < n and sample_intervals[j][1] <= c_start:
            j += 1
        k = j
        while k < n and sample_intervals[k][0] < c_end:
            s_start, s_end = sample_intervals[k]
            if s_end > c_start and s_start < c_end:
                vector[i] = 1
                break
            k += 1
    return vector


def build_support_matrix(
    samples: Sequence[Sample], consensus: Sequence[ConsensusInterval]
) -> Tuple[List[str], List[List[int]], List[int]]:
    """Return sample IDs, per-sample support vectors, and support counts."""
    sample_ids = [s.sample_id for s in samples]
    chr_to_consensus_indices: Dict[str, List[int]] = {}
    chr_to_consensus_intervals: Dict[str, List[Interval]] = {}
    for idx, (chrom, start, end) in enumerate(consensus):
        chr_to_consensus_indices.setdefault(chrom, []).append(idx)
        chr_to_consensus_intervals.setdefault(chrom, []).append((start, end))

    vectors: List[List[int]] = []
    for sample in samples:
        peaks = merge_by_chrom(read_bed3_like(sample.peak_path))
        sample_vec = [0] * len(consensus)
        for chrom, indices in chr_to_consensus_indices.items():
            c_intervals = chr_to_consensus_intervals[chrom]
            p_intervals = peaks.get(chrom, [])
            chrom_vec = intervals_support_vector(c_intervals, p_intervals)
            for local_i, global_i in enumerate(indices):
                sample_vec[global_i] = chrom_vec[local_i]
        vectors.append(sample_vec)

    support_counts = [sum(v[i] for v in vectors) for i in range(len(consensus))]
    return sample_ids, vectors, support_counts


def write_consensus_bed(consensus: Iterable[ConsensusInterval], output_path: str) -> None:
    with open(output_path, "w", encoding="utf-8") as out:
        for chrom, start, end in consensus:
            out.write(f"{chrom}\t{start}\t{end}\n")


def write_matrix_tsv(
    output_path: str,
    consensus: Sequence[ConsensusInterval],
    sample_ids: Sequence[str],
    sample_vectors: Sequence[Sequence[int]],
    support_counts: Sequence[int],
    header: bool = True,
) -> None:
    with open(output_path, "w", encoding="utf-8") as out:
        if header:
            out.write("chr\tstart\tend\tsupport_count\t" + "\t".join(sample_ids) + "\n")
        for i, (chrom, start, end) in enumerate(consensus):
            flags = [str(sample_vectors[s_idx][i]) for s_idx in range(len(sample_ids))]
            out.write(
                f"{chrom}\t{start}\t{end}\t{support_counts[i]}\t" + "\t".join(flags) + "\n"
            )
