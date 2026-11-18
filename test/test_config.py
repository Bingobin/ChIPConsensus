"""Regression tests for named sample config columns."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src.consensus_lib import Sample, load_samples_from_config


class LoadSamplesFromConfigTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.root = Path(self.temp_dir.name)
        self.peak = self.root / "sample, peaks.bed"
        self.peak.write_text("chr1\t10\t20\n", encoding="utf-8")
        self.config = self.root / "samples.csv"

    def write_config(self, rows, delimiter=",", encoding="utf-8", preamble=""):
        with self.config.open("w", encoding=encoding, newline="") as fh:
            fh.write(preamble)
            csv.writer(fh, delimiter=delimiter).writerows(rows)

    def assert_sample_loaded(self):
        self.assertEqual(
            load_samples_from_config(str(self.config)),
            [Sample(sample_id="sample1", peak_path=str(self.peak))],
        )

    def test_original_two_column_csv_and_tsv(self):
        for delimiter in (",", "\t"):
            with self.subTest(delimiter=delimiter):
                self.write_config([["ID", "PeakFile"], ["sample1", str(self.peak)]], delimiter)
                self.assert_sample_loaded()

    def test_extra_columns_are_ignored(self):
        for delimiter in (",", "\t"):
            with self.subTest(delimiter=delimiter):
                self.write_config(
                    [["ID", "PeakFile", "Big"], ["sample1", str(self.peak), "missing.bw"]],
                    delimiter,
                )
                self.assert_sample_loaded()

    def test_reordered_columns(self):
        for delimiter in (",", "\t"):
            with self.subTest(delimiter=delimiter):
                self.write_config(
                    [["Big", "PeakFile", "ID"], ["missing.bw", str(self.peak), "sample1"]],
                    delimiter,
                )
                self.assert_sample_loaded()

    def test_case_whitespace_bom_and_comments(self):
        self.write_config(
            [[" big ", " peakFILE ", " id "], ["", f" {self.peak} ", " sample1 "]],
            encoding="utf-8-sig",
            preamble="\n  # comment before header\n",
        )
        self.assert_sample_loaded()

    def test_omitted_optional_values(self):
        for delimiter in (",", "\t"):
            with self.subTest(delimiter=delimiter):
                self.write_config(
                    [["ID", "PeakFile", "Big"], ["sample1", str(self.peak)]], delimiter
                )
                self.assert_sample_loaded()

    def test_missing_required_header(self):
        cases = ((["ID", "Big"], "PeakFile"), (["Big", "PeakFile"], "ID"))
        for header, missing in cases:
            with self.subTest(missing=missing):
                self.write_config([header])
                with self.assertRaisesRegex(ValueError, f"missing required columns: {missing}"):
                    load_samples_from_config(str(self.config))

    def test_duplicate_required_header(self):
        self.write_config([["ID", "PeakFile", "id"]])
        with self.assertRaisesRegex(ValueError, "duplicate required column: ID"):
            load_samples_from_config(str(self.config))

    def test_empty_required_values(self):
        cases = ((["", str(self.peak), "x"], "ID"), (["sample1", "", "x"], "PeakFile"))
        for delimiter in (",", "\t"):
            for row, missing in cases:
                with self.subTest(delimiter=delimiter, missing=missing):
                    self.write_config([["ID", "PeakFile", "Big"], row], delimiter)
                    with self.assertRaisesRegex(ValueError, f"line 2: missing value for {missing}"):
                        load_samples_from_config(str(self.config))

    def test_missing_required_value_after_reordering(self):
        self.write_config([["Big", "PeakFile", "ID"], ["x", str(self.peak)]])
        with self.assertRaisesRegex(ValueError, "line 2: missing value for ID"):
            load_samples_from_config(str(self.config))

    def test_duplicate_sample_id(self):
        self.write_config(
            [["ID", "PeakFile"], ["sample1", str(self.peak)], ["sample1", str(self.peak)]]
        )
        with self.assertRaisesRegex(ValueError, "Duplicate sample ID 'sample1' at line 3"):
            load_samples_from_config(str(self.config))

    def test_missing_peak_file(self):
        self.write_config([["ID", "PeakFile"], ["sample1", str(self.root / "missing.bed")]])
        with self.assertRaisesRegex(FileNotFoundError, "Peak file not found for sample 'sample1'"):
            load_samples_from_config(str(self.config))

    def test_empty_config_and_no_samples(self):
        cases = (([], "No valid header found"), ([["ID", "PeakFile"]], "No samples found"))
        for rows, message in cases:
            with self.subTest(rows=rows):
                self.write_config(rows, preamble="# comment\n\n")
                with self.assertRaisesRegex(ValueError, message):
                    load_samples_from_config(str(self.config))

    def test_pipeline_accepts_extra_and_reordered_columns(self):
        repo = Path(__file__).resolve().parents[1]
        configs = (
            [["ID", "PeakFile", "Big"], ["sample1", str(self.peak), "missing.bw"]],
            [["Big", "PeakFile", "ID"], ["missing.bw", str(self.peak), "sample1"]],
        )
        for rows in configs:
            with self.subTest(header=rows[0]):
                self.write_config(rows)
                result = subprocess.run(
                    [sys.executable, "-B", str(repo / "chip_consensus.py"),
                     "--config", str(self.config), "--prefix", "test",
                     "--out-dir", str(self.root), "--no-upset"],
                    cwd=self.root, capture_output=True, text=True,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(
                    (self.root / "test.consensus.support.tsv").read_text(encoding="utf-8"),
                    "chr\tstart\tend\tsupport_count\tsample1\nchr1\t10\t20\t1\t1\n",
                )


if __name__ == "__main__":
    unittest.main()
