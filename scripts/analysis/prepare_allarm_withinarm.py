"""Prepare the primary all-arm treatment-by-weight dataset."""

from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[2]
SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-6.xlsx"
OUTPUTS = PACKAGE / "data" / "outputs"
DAYS_PER_MONTH = 30.4


def main():
    mice = pd.read_excel(SUPPLEMENT, sheet_name="SD1_mouse_level")
    meta = pd.read_csv(OUTPUTS / "withinarm_interactions.csv")
    for age in (6, 12, 18, 24):
        col = f"bw_{age}"
        mice[col] = pd.to_numeric(mice[col], errors="coerce")
        mice.loc[~mice[col].between(10, 80), col] = np.nan
    mice = mice[mice.dead.notna() & mice.lifespan_days.notna()].copy()

    parts = []
    for row in meta[["cohort", "group", "sex", "weight_age"]].drop_duplicates().itertuples(index=False):
        age = int(row.weight_age)
        weight_col = f"bw_{age}"
        entry = age * DAYS_PER_MONTH
        z = mice[
            (mice.cohort == row.cohort)
            & (mice.sex == row.sex)
            & mice.group.isin(["Control", row.group])
            & mice[weight_col].notna()
            & (mice.lifespan_days > entry)
        ].copy()
        z["Tx"] = (z.group == row.group).astype(int)
        if z.Tx.nunique() < 2 or len(z) < 60 or z.dead.sum() < 40:
            continue
        z["wz"] = z.groupby(["site", "Tx"])[weight_col].transform(
            lambda x: (x - x.mean()) / x.std(ddof=1)
        )
        z = z.dropna(subset=["wz"])
        z["entry"] = entry
        z["armsex"] = f"{row.cohort}|{row.group}|{row.sex}"
        z["stratum"] = z.site.astype(str) + "|" + z.armsex
        z["mouse"] = z.cohort.astype(str) + "|" + z.id.astype(str)
        z["interaction"] = z.Tx * z.wz
        parts.append(z[["sex", "armsex", "stratum", "mouse", "Tx", "wz", "entry", "lifespan_days", "dead", "interaction"]])

    out = pd.concat(parts, ignore_index=True)
    path = OUTPUTS / "itp_allarm_withinarm_long_reproduced.csv"
    out.to_csv(path, index=False)
    summary = out.groupby("sex").agg(arms=("armsex", "nunique"), rows=("armsex", "size"), unique_mice=("mouse", "nunique"))
    expected = {"m": (89, 14301), "f": (85, 12416)}
    for sex, (arms, unique_mice) in expected.items():
        got = summary.loc[sex]
        if int(got.arms) != arms or int(got.unique_mice) != unique_mice:
            raise RuntimeError(f"Primary reconciliation failed for {sex}: {got.to_dict()}")
    print(summary)
    print(path)


if __name__ == "__main__":
    main()
