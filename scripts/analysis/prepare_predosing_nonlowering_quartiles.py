"""Build the unique-mouse 29-arm pre-dosing non-lowering sensitivity set."""
from pathlib import Path

import pandas as pd


PACKAGE = Path(__file__).resolve().parents[2]
SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-6.xlsx"
INPUTS = PACKAGE / "data" / "inputs"
OUTPUTS = PACKAGE / "data" / "outputs"


def main() -> None:
    mice = pd.read_excel(SUPPLEMENT, sheet_name="SD1_mouse_level")
    weight_effects = pd.read_csv(OUTPUTS / "weight_effect_by_arm.csv")
    arm_meta = pd.read_csv(OUTPUTS / "withinarm_interactions.csv")

    post = weight_effects[
        weight_effects["post_dosing"].astype(str).str.lower().eq("true")
    ].copy()
    post["weighted_diff"] = post["diff_g"] * post["n_tx"]
    classification = post.groupby(["cohort", "group"], as_index=False).agg(
        weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum")
    )
    classification["mean_diff_g"] = (
        classification["weighted_diff"] / classification["n_tx"]
    )
    initiation = post[["cohort", "group", "init"]].drop_duplicates()
    eligible = arm_meta.merge(
        initiation, on=["cohort", "group"], how="left"
    ).merge(
        classification[["cohort", "group", "mean_diff_g"]],
        on=["cohort", "group"], how="left",
    )
    eligible = eligible[
        (eligible["sex"] == "m")
        & (eligible["mean_diff_g"] >= -1)
        & (eligible["weight_age"] < eligible["init"])
    ].drop_duplicates(["cohort", "group"])

    parts = []
    for row in eligible.itertuples(index=False):
        weight_col = f"bw_{int(row.weight_age)}"
        entry = float(row.weight_age) * 30.4375
        subset = mice[
            (mice["cohort"] == row.cohort)
            & (mice["sex"] == "m")
            & mice["group"].isin(["Control", row.group])
            & mice[weight_col].notna()
            & (mice["lifespan_days"] > entry)
        ].copy()
        subset["Tx"] = (subset["group"] == row.group).astype(int)
        subset["entry"] = entry
        subset["body_weight"] = pd.to_numeric(subset[weight_col], errors="coerce")
        parts.append(subset)

    data = pd.concat(parts, ignore_index=True)
    data["mouse"] = data["cohort"].astype(str) + "|" + data["id"].astype(str)
    data = data.drop_duplicates("mouse")
    if len(eligible) != 29 or len(data) != 6658:
        raise RuntimeError(
            f"Expected 29 arms/6,658 mice; found {len(eligible)}/{len(data)}"
        )

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    output = OUTPUTS / "predosing_weight_nonlowering_29arm_unique.csv"
    columns = [
        "cohort", "site", "id", "mouse", "group", "Tx", "entry",
        "body_weight", "lifespan_days", "dead",
    ]
    data[columns].to_csv(output, index=False)
    print({"arms": len(eligible), "unique_mice": len(data), "output": str(output)})


if __name__ == "__main__":
    main()
