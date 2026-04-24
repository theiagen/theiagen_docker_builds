#!/usr/bin/env python3
"""Transform freyja variant abundances to long format with per-site/date averaged variant percents."""

import argparse
import sys
import pandas as pd

REQUIRED_META_COLS = ["collection_site", "collection_date", "month", "day", "year"]
OPTIONAL_META_COLS = ["latitude", "longitude"]
METADATA_COLS = REQUIRED_META_COLS + OPTIONAL_META_COLS
GROUP_KEY_COLS = ["collection_site", "collection_date"]  # uniquely identifies a site+date


def normalize_percents(x):
    """Round percents in a group so they sum to exactly 100."""
    x = x.round(4).copy()
    x.iloc[-1] = round(100 - x.iloc[:-1].sum(), 4)
    return x


def derive_date_parts(df):
    """Parse month, day, year from collection_date if not already present."""
    parsed = pd.to_datetime(df["collection_date"])
    for col, attr in [("month", "month"), ("day", "day"), ("year", "year")]:
        if col not in df.columns:
            df[col] = getattr(parsed.dt, attr)
    return df


def compute_epiweek(dates):
    """MMWR epiweek label (Sunday-Saturday) formatted as 'YYYYWww'."""
    # Shift Sunday-start weeks to Monday-start so ISO week numbering gives MMWR weeks.
    shifted = dates + pd.Timedelta(days=1)
    iso = shifted.dt.isocalendar()
    return iso["year"].astype(str) + "W" + iso["week"].astype(str).str.zfill(2)


def freyja_to_long(input_tsv, output_csv, sample_col, group_by="date"):
    required_cols = [sample_col, "freyja_lineages", "freyja_abundances"] + METADATA_COLS
    existing_cols = set(pd.read_csv(input_tsv, sep="\t", nrows=0).columns)
    freyja_variance_df = pd.read_csv(
        input_tsv,
        sep="\t",
        usecols=[c for c in required_cols if c in existing_cols],
        dtype={"freyja_lineages": str, "freyja_abundances": str},
    )

    if group_by == "week":
        freyja_variance_df["week"] = compute_epiweek(pd.to_datetime(freyja_variance_df["collection_date"]))
        group_key_cols = ["collection_site", "week"]
        meta_cols = group_key_cols + [c for c in OPTIONAL_META_COLS if c in existing_cols]
    else:
        freyja_variance_df = derive_date_parts(freyja_variance_df)
        group_key_cols = GROUP_KEY_COLS
        meta_cols = [c for c in METADATA_COLS if c in freyja_variance_df.columns]

    # drop rows missing lineage/abundance data
    freyja_variance_clean_df = freyja_variance_df.dropna(subset=["freyja_lineages", "freyja_abundances"]).copy()

    # warn and drop rows where lineage/abundance counts don't match
    lineage_counts = freyja_variance_clean_df["freyja_lineages"].str.split(",").str.len()
    abundance_counts = freyja_variance_clean_df["freyja_abundances"].str.split(",").str.len()
    mismatch = lineage_counts != abundance_counts
    if mismatch.any():
        print(f"Warning: skipping {mismatch.sum()} rows with lineage/abundance count mismatch", file=sys.stderr)
        freyja_variance_clean_df = freyja_variance_clean_df[~mismatch]

    # split comma-separated strings into lists, then explode to long format
    freyja_variance_clean_df["variant"] = freyja_variance_clean_df["freyja_lineages"].str.split(",")
    freyja_variance_clean_df["percent"] = freyja_variance_clean_df["freyja_abundances"].str.split(",")
    freyja_variance_long_format_df = freyja_variance_clean_df.explode(["variant", "percent"]).reset_index(drop=True)
    freyja_variance_long_format_df["variant"] = freyja_variance_long_format_df["variant"].str.strip()
    freyja_variance_long_format_df["percent"] = freyja_variance_long_format_df["percent"].astype(float)

    # normalize per sample so each sample's variants sum to 100%
    freyja_variance_long_format_df["percent"] = freyja_variance_long_format_df.groupby(sample_col)["percent"].transform(lambda x: x / x.sum() * 100)

    # count samples per group before aggregating (absent variants contribute 0)
    sample_counts = freyja_variance_long_format_df.groupby(group_key_cols)[sample_col].nunique().rename("n_samples")

    # sum per group+variant, divide by total samples in group; collect contributing samples
    present_meta = [c for c in meta_cols if c in freyja_variance_long_format_df.columns]
    extra_meta = [c for c in present_meta if c not in group_key_cols]
    freyja_variance_long_format_df = freyja_variance_long_format_df.groupby(group_key_cols + extra_meta + ["variant"], as_index=False).agg(
        percent=("percent", "sum"),
        origin_samples=(sample_col, lambda x: ",".join(sorted(x.unique()))),
    )
    freyja_variance_long_format_df = freyja_variance_long_format_df.merge(sample_counts, on=group_key_cols)
    freyja_variance_long_format_df["percent"] = freyja_variance_long_format_df["percent"] / freyja_variance_long_format_df.pop("n_samples")

    # fix rounding so each group sums to exactly 100
    freyja_variance_long_format_df["percent"] = freyja_variance_long_format_df.groupby(group_key_cols)["percent"].transform(normalize_percents)

    output_cols = ["id"] + present_meta[:1] + ["variant", "percent", "origin_samples"] + present_meta[1:]
    freyja_variance_long_format_df.insert(0, "id", range(1, len(freyja_variance_long_format_df) + 1))
    freyja_variance_long_format_df[output_cols].to_csv(output_csv, sep="\t", index=False)
    print(f"Wrote {len(freyja_variance_long_format_df)} rows to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transform Freyja TSV to long-format TSV with averaged variant percents per site/date.")
    parser.add_argument("input_tsv", help="Input Freyja TSV file")
    parser.add_argument("output_tsv", help="Output long-format TSV file")
    parser.add_argument("--sample-col", required=True, help="Column name in the TSV to use as the sample identifier")
    parser.add_argument("--group-by", choices=["date", "week"], default="date",
                        help="Aggregation key: 'date' groups by site+collection_date (default); 'week' groups by site+MMWR epiweek (Sun-Sat) and drops date parts from output")
    args = parser.parse_args()
    freyja_to_long(args.input_tsv, args.output_tsv, args.sample_col, args.group_by)
