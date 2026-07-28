#!/usr/bin/env python3
"""
merge_anubis_tables.py — merge ANUBIS result tables from multiple samples.

Usage
-----
python merge_anubis_tables.py sample1/anubis_results.tsv sample2/anubis_results.tsv \
    -o merged_anubis.tsv

Each input file is expected to be an `anubis_results.tsv` produced by `anubis`.
A `sample` column is added (inferred from the parent directory name) and all
tables are concatenated into one long-format TSV.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    sample = path.parent.name if path.parent.name else path.stem
    df.insert(0, "sample", sample)
    return df


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tables", nargs="+", type=Path,
                   help="anubis_results.tsv files (one per sample)")
    p.add_argument("-o", "--output", default="merged_anubis.tsv",
                   help="Output TSV path (default: merged_anubis.tsv)")
    p.add_argument("--sample-names", nargs="+",
                   help="Override sample names (must match number of input files)")
    args = p.parse_args(argv)

    if args.sample_names and len(args.sample_names) != len(args.tables):
        print("ERROR: --sample-names count must match number of input tables",
              file=sys.stderr)
        sys.exit(1)

    frames = []
    for i, path in enumerate(args.tables):
        df = load_table(path)
        if args.sample_names:
            df["sample"] = args.sample_names[i]
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Merged {len(frames)} tables ({len(merged)} rows) → {args.output}")


if __name__ == "__main__":
    main()
