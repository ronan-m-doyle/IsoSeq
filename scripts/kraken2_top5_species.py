#!/usr/bin/env python3

import pandas as pd
import glob
import os
import re
import argparse
import sys


def parse_kraken_report(path):
    rows = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 6:
                percent, reads_clade, reads_direct, rank, taxid = parts[:5]
                name = '\t'.join(parts[5:]).strip()
            else:
                parts_ws = re.split(r'\s+', line, maxsplit=5)
                if len(parts_ws) < 6:
                    continue
                percent, reads_clade, reads_direct, rank, taxid, name = parts_ws
                name = name.strip()

            try:
                percent_f = float(percent.strip().rstrip('%'))
            except Exception:
                continue

            try:
                reads_clade_i = int(float(reads_clade.strip()))
            except Exception:
                reads_clade_i = 0

            try:
                reads_direct_i = int(float(reads_direct.strip()))
            except Exception:
                reads_direct_i = 0

            rows.append({
                "percent": percent_f,
                "reads_clade": reads_clade_i,
                "reads_direct": reads_direct_i,
                "rank": rank.strip(),
                "taxid": taxid.strip(),
                "name": name
            })

    return pd.DataFrame(rows)


# ---------------- CLI ----------------

parser = argparse.ArgumentParser(
    description="Extract top 5 species from Kraken2 report(s)."
)

parser.add_argument(
    "--input",
    required=True,
    help="Input kraken2_report.tsv file or glob pattern (e.g. '*/kraken2_report.tsv')"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output CSV file name"
)

args = parser.parse_args()

files = glob.glob(args.input)

if not files:
    print(f"❌ No files found matching: {args.input}")
    sys.exit(1)

all_top5 = []

for path in sorted(files):
    sample = os.path.basename(os.path.dirname(path))
    df = parse_kraken_report(path)

    if df.empty:
        continue

    df_species = df[df["rank"] == "S"].copy()
    if df_species.empty:
        continue

    df_species["percent"] = pd.to_numeric(df_species["percent"], errors="coerce").fillna(0.0)
    df_species["reads_clade"] = pd.to_numeric(df_species["reads_clade"], errors="coerce").fillna(0).astype(int)

    df_top5 = (
        df_species
        .sort_values("percent", ascending=False)
        .head(5)
        [["name", "percent", "reads_clade", "taxid"]]
        .copy()
    )

    df_top5["sample"] = sample
    df_top5 = df_top5[["sample", "name", "percent", "reads_clade", "taxid"]]

    all_top5.append(df_top5)

if all_top5:
    summary = (
        pd.concat(all_top5, ignore_index=True)
        .sort_values(["sample", "percent"], ascending=[True, False])
        .reset_index(drop=True)
    )

    summary.to_csv(args.output, index=False)
    print(f"✅ Saved {args.output}")
else:
    print("No valid species-level entries found.")
    sys.exit(1)