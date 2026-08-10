"""Reproduce Table 4's descriptive comparator-state tertiles.

Each Kaplan-Meier median is left-truncated at six months. Percentage gain is
calculated within every cohort x site x sex x intervention cell, then averaged
without weights within sex-specific comparator-weight tertiles.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter


ROOT = Path(__file__).resolve().parents[2]
LIFE = ROOT / "data/inputs/ITP_bodyweight_all_cohorts.xlsx"
CELLS = ROOT / "data/outputs/comparator_effect_cells.csv"
SUMMARY_OUT = ROOT / "data/outputs/comparator_state_reproduced.csv"
CELL_OUT = ROOT / "data/outputs/comparator_cell_medians.csv"
Z95 = 1.959963984540054
LANDMARK_DAYS = 6 * 30.4


def km_median(frame):
    """Six-month-left-truncated Kaplan-Meier median."""
    z = frame.loc[frame["lifespan_days"].gt(LANDMARK_DAYS),
                  ["lifespan_days", "dead"]].dropna()
    if z.empty:
        return np.nan
    fit = KaplanMeierFitter()
    fit.fit(
        z["lifespan_days"],
        event_observed=z["dead"],
        entry=np.full(len(z), LANDMARK_DAYS),
    )
    return float(fit.median_survival_time_)


life = pd.read_excel(LIFE)
for column in ("lifespan_days", "dead"):
    life[column] = pd.to_numeric(life[column], errors="coerce")
life = life[life.lifespan_days.notna() & life.dead.notna()].copy()
cells = pd.read_csv(CELLS)
if len(cells) != 426:
    raise AssertionError(f"Expected 426 comparator cells, found {len(cells)}")

medians = []
for row in cells.itertuples(index=False):
    treated = life[(life.cohort == row.cohort) & (life.group == row.group) &
                   (life.site == row.site) & (life.sex == row.sex)]
    control = life[(life.cohort == row.cohort) & (life.group == "Control") &
                   (life.site == row.site) & (life.sex == row.sex)]
    cmed, tmed = km_median(control), km_median(treated)
    medians.append((cmed, tmed, 100 * (tmed - cmed) / cmed))

cells[["control_km_median", "treated_km_median", "median_gain_pct"]] = pd.DataFrame(
    medians, index=cells.index
)
cells["control_cluster"] = cells.cohort.astype(str) + "|" + cells.site.astype(str)
cells["sex_specific_control"] = cells.control_cluster + "|" + cells.sex.astype(str)
cells.to_csv(CELL_OUT, index=False)

rows = []
for sex in ("m", "f"):
    subset = cells[cells.sex == sex].copy()
    for split in (3, 4):
        # Match the authoritative analysis: pandas qcut ranks the observed
        # control-weight values separately within sex.
        subset["bin"] = pd.qcut(
            subset.ctl_bw6, split, labels=[f"Q{i}" for i in range(1, split + 1)]
        )
        for label, group in subset.groupby("bin", observed=True):
            weights = 1 / group.se.pow(2)
            beta = np.average(group.log_hr, weights=weights)
            se = np.sqrt(1 / weights.sum())
            rows.append({
                "sex": sex,
                "split": split,
                "bin": str(label),
                "cells": len(group),
                "mean_control_weight_g": group.ctl_bw6.mean(),
                "pooled_hr": np.exp(beta),
                "lower_95_ci": np.exp(beta - Z95 * se),
                "upper_95_ci": np.exp(beta + Z95 * se),
                "hazard_reduction_pct": 100 * (1 - np.exp(beta)),
                "mean_cell_km_median_gain_pct": group.median_gain_pct.mean(),
            })

summary = pd.DataFrame(rows)
summary.to_csv(SUMMARY_OUT, index=False)

print(f"{len(cells)} cells; {cells.control_cluster.nunique()} cohort-site control clusters; "
      f"{cells.sex_specific_control.nunique()} sex-specific control means")
print(summary.to_string(index=False))

t3 = summary[summary.split == 3].set_index(["sex", "bin"])
expected = {
    ("m", "Q1", "hazard_reduction_pct"): 11.0,
    ("m", "Q3", "hazard_reduction_pct"): 22.3,
    ("f", "Q1", "hazard_reduction_pct"): 4.3,
    ("f", "Q3", "hazard_reduction_pct"): 15.1,
    ("m", "Q1", "mean_cell_km_median_gain_pct"): 4.1,
    ("m", "Q3", "mean_cell_km_median_gain_pct"): 9.3,
    ("f", "Q1", "mean_cell_km_median_gain_pct"): 1.6,
    ("f", "Q3", "mean_cell_km_median_gain_pct"): 4.1,
}
for (sex, bin_, field), wanted in expected.items():
    got = round(float(t3.loc[(sex, bin_), field]), 1)
    if got != wanted:
        raise AssertionError(f"{sex} {bin_} {field}: expected {wanted}, got {got}")
