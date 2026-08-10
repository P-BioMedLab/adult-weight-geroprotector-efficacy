"""Derive arm-level analysis tables from the harmonized mouse-level file."""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data" / "inputs" / "ITP_bodyweight_all_cohorts.xlsx"
OUTPUT = ROOT / "data" / "outputs"
DAYS_PER_MONTH = 30.4
WEIGHTS = {6: "bw_6", 12: "bw_12", 18: "bw_18", 24: "bw_24"}


def standardize(values: pd.Series) -> pd.Series:
    return (values - values.mean()) / values.std(ddof=1)


def load_mouse_data() -> pd.DataFrame:
    mice = pd.read_excel(SOURCE)
    for column in WEIGHTS.values():
        mice[column] = pd.to_numeric(mice[column], errors="coerce")
        mice.loc[(mice[column] < 10) | (mice[column] > 80), column] = np.nan
    mice["init"] = pd.to_numeric(
        mice["age_initiation"].astype(str).str.split().str[0], errors="coerce"
    )
    return mice[mice["dead"].notna() & mice["lifespan_days"].notna()].copy()


def derive_weight_effects(mice: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    treated = mice[mice["group"] != "Control"]
    for (cohort, group), arm in treated.groupby(["cohort", "group"]):
        initiation = arm["init"].iloc[0]
        controls = mice[(mice["cohort"] == cohort) & (mice["group"] == "Control")]
        for age, column in WEIGHTS.items():
            tx_all = arm[column].dropna()
            ctl_all = controls[column].dropna()
            if len(tx_all) < 40 or len(ctl_all) < 40:
                continue
            differences, sample_sizes = [], []
            for sex in ["m", "f"]:
                tx = arm.loc[arm["sex"] == sex, column].dropna()
                ctl = controls.loc[controls["sex"] == sex, column].dropna()
                if len(tx) < 20 or len(ctl) < 20:
                    continue
                differences.append(tx.mean() - ctl.mean())
                sample_sizes.append(len(tx))
            if not differences:
                continue
            difference = float(np.average(differences, weights=sample_sizes))
            rows.append({
                "cohort": cohort,
                "group": group,
                "init": initiation,
                "weight_age": age,
                "post_dosing": bool(not np.isnan(initiation) and age >= initiation),
                "n_tx": len(tx_all),
                "diff_g": difference,
                "pct": 100 * difference / ctl_all.mean(),
                "p": ttest_ind(tx_all, ctl_all, equal_var=False).pvalue,
            })
    return pd.DataFrame(rows)


def fit_arm_interaction(
    mice: pd.DataFrame, cohort: str, group: str, age: int
) -> list[dict]:
    column = WEIGHTS[age]
    landmark = age * DAYS_PER_MONTH
    arm = mice[(mice["cohort"] == cohort) & mice["group"].isin([group, "Control"])].copy()
    arm = arm[arm[column].notna() & (arm["lifespan_days"] > landmark)]
    arm["Tx"] = (arm["group"] == group).astype(int)
    rows: list[dict] = []
    for sex in ["m", "f"]:
        frame = arm[arm["sex"] == sex].copy()
        if frame["Tx"].nunique() < 2 or len(frame) < 60 or frame["dead"].sum() < 40:
            continue
        frame["wz"] = frame.groupby(["site", "Tx"])[column].transform(standardize)
        frame = frame.dropna(subset=["wz"])
        frame["entry"] = landmark
        frame["Tx_x_wz"] = frame["Tx"] * frame["wz"]
        try:
            model = CoxPHFitter().fit(
                frame[["lifespan_days", "dead", "entry", "Tx", "wz", "Tx_x_wz", "site"]],
                "lifespan_days", "dead", entry_col="entry", strata=["site"],
            )
        except Exception:
            continue
        result = model.summary.loc["Tx_x_wz"]
        rows.append({
            "cohort": cohort,
            "group": group,
            "sex": sex,
            "weight_age": age,
            "within_arm": True,
            "n": len(frame),
            "int_hr": result["exp(coef)"],
            "int_lo": result["exp(coef) lower 95%"],
            "int_hi": result["exp(coef) upper 95%"],
            "int_p": result["p"],
            "coef": result["coef"],
            "se": result["se(coef)"],
        })
    return rows


def derive_withinarm_interactions(mice: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    arms = mice[mice["group"] != "Control"][["cohort", "group"]].drop_duplicates()
    for arm in arms.itertuples(index=False):
        initiation = mice.loc[
            (mice["cohort"] == arm.cohort) & (mice["group"] == arm.group), "init"
        ].iloc[0]
        age = 6 if np.isnan(initiation) or initiation > 6 else 12
        rows.extend(fit_arm_interaction(mice, arm.cohort, arm.group, age))
    return pd.DataFrame(rows)


def derive_comparator_cells(mice: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    treated = mice[mice["group"] != "Control"]
    landmark = 6 * DAYS_PER_MONTH
    for (cohort, group), arm in treated.groupby(["cohort", "group"]):
        for site in arm["site"].dropna().unique():
            for sex in ["m", "f"]:
                tx = arm[(arm["site"] == site) & (arm["sex"] == sex)]
                ctl = mice[
                    (mice["cohort"] == cohort)
                    & (mice["group"] == "Control")
                    & (mice["site"] == site)
                    & (mice["sex"] == sex)
                ]
                if len(tx) < 25 or len(ctl) < 25:
                    continue
                frame = pd.concat([tx.assign(Tx=1), ctl.assign(Tx=0)])
                frame = frame[frame["lifespan_days"] > landmark].copy()
                if frame["dead"].sum() < 30:
                    continue
                frame["entry"] = landmark
                try:
                    model = CoxPHFitter().fit(
                        frame[["lifespan_days", "dead", "entry", "Tx"]],
                        "lifespan_days", "dead", entry_col="entry",
                    )
                except Exception:
                    continue
                result = model.summary.loc["Tx"]
                rows.append({
                    "cohort": cohort,
                    "group": group,
                    "site": site,
                    "sex": sex,
                    "ctl_bw6": ctl["bw_6"].mean(),
                    "ctl_median": ctl["lifespan_days"].median(),
                    "n_ctl": len(ctl),
                    "n_tx": len(tx),
                    "log_hr": result["coef"],
                    "se": result["se(coef)"],
                    "hr": result["exp(coef)"],
                })
    return pd.DataFrame(rows).dropna(subset=["ctl_bw6", "ctl_median"])


def main() -> None:
    warnings.filterwarnings("ignore")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mice = load_mouse_data()
    tables = {
        "weight_effect_by_arm.csv": derive_weight_effects(mice),
        "withinarm_interactions.csv": derive_withinarm_interactions(mice),
        "comparator_effect_cells.csv": derive_comparator_cells(mice),
    }
    for name, table in tables.items():
        table.to_csv(OUTPUT / name, index=False)
        print(f"{name}: {len(table):,} rows")


if __name__ == "__main__":
    main()
