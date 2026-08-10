"""Reproduce stratum-level consistency of the six-month weight-lifespan relation.

Within each cohort x intervention x sex stratum, site-specific Spearman
correlations are Fisher-z transformed and pooled with n-3 weights. Only natural
deaths with a recorded six-month weight are used, and site cells require at
least 12 animals. The output supports the concise descriptive robustness result
reported in manuscript v12.
"""

from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[2]
INPUT = PACKAGE / "supplementary" / "Supplementary_Data_1-6.xlsx"
OUTPUT = PACKAGE / "data" / "outputs" / "rederived_bw6_strata.csv"


def main() -> None:
    data = pd.read_excel(INPUT, sheet_name="SD1_mouse_level")
    data = data[
        data["dead"].eq(1)
        & data["bw_6"].notna()
        & data["lifespan_days"].notna()
    ].copy()

    rows = []
    for (cohort, sex, group), stratum in data.groupby(
        ["cohort", "sex", "group"], sort=True
    ):
        site_estimates = []
        for _, cell in stratum.groupby("site", sort=True):
            if len(cell) < 12:
                continue
            rho = cell[["bw_6", "lifespan_days"]].corr(method="spearman").iloc[0, 1]
            if not np.isfinite(rho):
                continue
            site_estimates.append(
                (np.arctanh(np.clip(rho, -0.999999, 0.999999)), len(cell) - 3, len(cell))
            )
        if not site_estimates:
            continue
        z = np.asarray([value[0] for value in site_estimates])
        weights = np.asarray([value[1] for value in site_estimates], dtype=float)
        rows.append(
            {
                "cohort": cohort,
                "sex": sex,
                "group": group,
                "nsite": len(site_estimates),
                "n": int(sum(value[2] for value in site_estimates)),
                "rho": float(np.tanh(np.average(z, weights=weights))),
            }
        )

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT, index=False)

    expected = {"m": (105, 104), "f": (101, 83)}
    for sex, (n_expected, negative_expected) in expected.items():
        subset = result[result["sex"].eq(sex)]
        observed = (len(subset), int(subset["rho"].lt(0).sum()))
        if observed != (n_expected, negative_expected):
            raise RuntimeError(
                f"Reconciliation failed for sex={sex}: {observed} != "
                f"{(n_expected, negative_expected)}"
            )
        print(
            f"{sex}: {observed[1]}/{observed[0]} negative; "
            f"median rho={subset['rho'].median():.3f}"
        )
    print(OUTPUT)


if __name__ == "__main__":
    main()
