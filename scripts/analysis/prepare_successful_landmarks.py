"""Prepare 6/12/18/24-month landmarks for published lifespan extenders."""

from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[2]
SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-6.xlsx"
INPUTS = PACKAGE / "data" / "inputs"
OUTPUTS = PACKAGE / "data" / "outputs"
DAYS = 30.4375


def main():
    mice = pd.read_excel(SUPPLEMENT, sheet_name="SD1_mouse_level")
    specs = pd.read_csv(INPUTS / "published_lifespan_extenders.csv")
    rows = []
    for spec in specs.itertuples(index=False):
        for landmark in (6, 12, 18, 24):
            if landmark <= spec.start_age:
                continue
            col = f"bw_{landmark}"
            z = mice[
                (mice.cohort == spec.cohort)
                & mice.group.isin(["Control", spec.group])
                & mice[col].notna()
                & (mice.lifespan_days > landmark * DAYS)
            ].copy()
            z["Tx"] = (z.group == spec.group).astype(int)
            z["bw"] = pd.to_numeric(z[col], errors="coerce")
            z = z[z.bw.between(10, 80)]
            z["weight_z"] = z.groupby(["site", "sex", "Tx"])["bw"].transform(
                lambda x: (x - x.mean()) / x.std(ddof=1)
            )
            z = z.dropna(subset=["weight_z"])
            z["entry"] = landmark * DAYS
            z["arm"] = f"{spec.cohort}|{spec.group}"
            z["arm_stage"] = z.arm + f"|L{landmark}"
            z["landmark"] = landmark
            z["positive_m"] = bool(spec.positive_m)
            z["positive_f"] = bool(spec.positive_f)
            z["mouse"] = z.cohort.astype(str) + "|" + z.id.astype(str)
            rows.append(z)
    out = pd.concat(rows, ignore_index=True)
    cols = ["cohort", "site", "sex", "id", "mouse", "group", "Tx", "dead", "lifespan_days", "entry", "bw", "weight_z", "arm", "arm_stage", "landmark", "positive_m", "positive_f"]
    path = OUTPUTS / "successful_landmarks_long.csv"
    out[cols].to_csv(path, index=False)
    print(f"{out.arm.nunique()} arms; {out.mouse.nunique()} unique mice; {len(out)} rows")
    print(path)


if __name__ == "__main__":
    main()
