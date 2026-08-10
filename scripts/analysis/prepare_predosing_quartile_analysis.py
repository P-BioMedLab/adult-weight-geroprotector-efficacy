"""Build the exact 27-arm pre-dosing weight-lowering comparison dataset.

Inputs are resolved relative to the packaged supplement and derived analysis
tables. Controls are expanded once per randomized intervention comparison;
mouse IDs allow the downstream Cox models to cluster reused controls.
"""

from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[2]
SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-6.xlsx"
INPUTS = PACKAGE / "data" / "inputs"
OUTPUTS = PACKAGE / "data" / "outputs"
DAYS_PER_MONTH = 30.4375


def quartile_from_rank(series: pd.Series) -> pd.Series:
    rank = series.rank(method="first", pct=True)
    return np.ceil(rank * 4).clip(1, 4).astype(int)


def main():
    mice = pd.read_excel(SUPPLEMENT, sheet_name="SD1_mouse_level")
    weight_effects = pd.read_csv(OUTPUTS / "weight_effect_by_arm.csv")
    arm_meta = pd.read_csv(OUTPUTS / "withinarm_interactions.csv")

    post = weight_effects[
        weight_effects["post_dosing"].astype(str).str.lower().eq("true")
    ].copy()
    post["weighted_diff"] = post["diff_g"] * post["n_tx"]
    classification = (
        post.groupby(["cohort", "group"], as_index=False)
        .agg(weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum"))
    )
    classification["mean_diff_g"] = (
        classification["weighted_diff"] / classification["n_tx"]
    )
    classification["weight_lowering"] = classification["mean_diff_g"] < -1

    init = post[["cohort", "group", "init"]].drop_duplicates()
    eligible = arm_meta.merge(init, on=["cohort", "group"], how="left")
    eligible = eligible.merge(
        classification[["cohort", "group", "mean_diff_g", "weight_lowering"]],
        on=["cohort", "group"],
        how="left",
    )
    eligible = eligible[
        (eligible["sex"] == "m")
        & eligible["weight_lowering"].eq(True)
        & (eligible["weight_age"] < eligible["init"])
    ].drop_duplicates(["cohort", "group"])

    parts = []
    for row in eligible.itertuples(index=False):
        weight_col = f"bw_{int(row.weight_age)}"
        entry = float(row.weight_age) * DAYS_PER_MONTH
        z = mice[
            (mice["cohort"] == row.cohort)
            & (mice["sex"] == "m")
            & mice["group"].isin(["Control", row.group])
            & mice[weight_col].notna()
            & (mice["lifespan_days"] > entry)
        ].copy()
        z["arm"] = f"{row.cohort}|{row.group}"
        z["Tx"] = (z["group"] == row.group).astype(int)
        z["weight_age"] = int(row.weight_age)
        z["entry"] = entry
        z["body_weight"] = pd.to_numeric(z[weight_col], errors="coerce")
        z["mean_treatment_weight_effect_g"] = row.mean_diff_g
        parts.append(z)

    data = pd.concat(parts, ignore_index=True)
    data["mouse"] = data["cohort"].astype(str) + "|" + data["id"].astype(str)
    data["stratum"] = data["arm"] + "|" + data["site"].astype(str)
    data["matched_cell"] = data["stratum"] + "|Tx" + data["Tx"].astype(str)
    data["matched_quartile"] = data.groupby("matched_cell", group_keys=False)[
        "body_weight"
    ].apply(quartile_from_rank)
    data["weight_z_matched"] = data.groupby("matched_cell")["body_weight"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=1)
    )

    # Direct standardization: every arm-site-treatment cell contributes the
    # same total weight within each quartile. This removes changes in arm/site
    # composition from the displayed survival curves and medians.
    n_cell_q = data.groupby(["matched_cell", "matched_quartile"])["mouse"].transform("size")
    data["standardization_weight"] = 1.0 / n_cell_q

    # Reused controls are identifiable and handled by mouse-clustered variance.
    data["comparison_reuse_count"] = data.groupby("mouse")["arm"].transform("nunique")

    cols = [
        "cohort", "site", "id", "mouse", "arm", "stratum", "group", "Tx",
        "weight_age", "entry", "body_weight", "weight_z_matched",
        "matched_quartile", "standardization_weight", "comparison_reuse_count",
        "lifespan_days", "dead", "mean_treatment_weight_effect_g",
    ]
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    out = OUTPUTS / "predosing_weight_lowering_27arm_comparisons.csv"
    data[cols].to_csv(out, index=False)
    eligible[["cohort", "group", "weight_age", "init", "mean_diff_g"]].to_csv(
        OUTPUTS / "predosing_weight_lowering_27arm_eligibility.csv", index=False
    )

    summary = {
        "arms": int(data["arm"].nunique()),
        "comparison_rows": int(len(data)),
        "unique_mice": int(data["mouse"].nunique()),
        "treated_mice": int(data.loc[data.Tx.eq(1), "mouse"].nunique()),
        "control_mice": int(data.loc[data.Tx.eq(0), "mouse"].nunique()),
    }
    expected = {"arms": 27, "comparison_rows": 10018, "unique_mice": 6080}
    for key, value in expected.items():
        if summary[key] != value:
            raise RuntimeError(f"Reconciliation failed for {key}: {summary[key]} != {value}")
    print(summary)
    print(out)


if __name__ == "__main__":
    main()
