"""Reproduce the cross-study and body-composition summaries.

Run from any directory. Inputs and outputs are resolved relative to this package.
Requires numpy, pandas and an Excel engine supported by pandas.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PACKAGE = Path(__file__).resolve().parents[2]
SOURCE = PACKAGE / "supplementary" / "Supplementary_Data_1-6.xlsx"
OUTPUT = PACKAGE / "data" / "outputs" / "supplementary_analyses_reproduced.csv"


def clustered_ols(X, y, cluster):
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    beta = np.linalg.solve(X.T @ X, X.T @ y)
    resid = y - X @ beta
    groups = np.unique(cluster)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for group in groups:
        score = X[cluster == group].T @ resid[cluster == group]
        meat += np.outer(score, score)
    bread = np.linalg.inv(X.T @ X)
    n, k, g = len(y), X.shape[1], len(groups)
    correction = (g / (g - 1)) * ((n - 1) / (n - k))
    return beta, correction * bread @ meat @ bread


def normal_p(z):
    return math.erfc(abs(float(z)) / math.sqrt(2))


Z = stats.norm.ppf(0.975)


def fisher_ci(r, n):
    transformed, se = math.atanh(r), 1 / math.sqrt(n - 3)
    return math.tanh(transformed - Z * se), math.tanh(transformed + Z * se)


rows = []

# SD3: responder analysis (weight decreased and lifespan increased).
d = pd.read_excel(SOURCE, sheet_name="SD3_cross_study_arms")
d = d[(d.weight_change_percent < 0) & (d.lifespan_change_percent > 0)]
d = d[d.gender.isin(["Male", "Female"])].copy()
d["x"] = -d.weight_change_percent.astype(float)
d["y"] = d.lifespan_change_percent.astype(float)
d["study"] = d.pubmed_id.astype(str)

for (source, sex), z in d.groupby(["source_group", "gender"]):
    beta, vcov = clustered_ols(z[["x"]], z.y, z.study.to_numpy())
    se = math.sqrt(vcov[0, 0])
    rows.append([f"slope_{source}_{sex}", len(z), z.study.nunique(), beta[0], beta[0] - Z * se, beta[0] + Z * se, normal_p(beta[0] / se)])
    rows.append([f"median_ratio_{source}_{sex}", len(z), z.study.nunique(),
                 (z.y / z.x).median(), np.nan, np.nan, np.nan])

for source in ["ITP", "non-ITP pharmacological", "dietary restriction"]:
    z = d[d.source_group == source].copy()
    female = (z.gender == "Female").astype(float).to_numpy()
    X = np.column_stack([z.x.to_numpy(), z.x.to_numpy() * female])
    beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
    se = math.sqrt(vcov[1, 1])
    rows.append([f"female_minus_male_{source}", len(z), z.study.nunique(), beta[1], beta[1] - Z * se, beta[1] + Z * se, normal_p(beta[1] / se)])

z = d[d.source_group.isin(["ITP", "non-ITP pharmacological"])].copy()
female = (z.gender == "Female").astype(float).to_numpy()
non_itp = (z.source_group == "non-ITP pharmacological").astype(float).to_numpy()
x = z.x.to_numpy()
X = np.column_stack([x, x * female, x * non_itp, x * female * non_itp])
beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
se = math.sqrt(vcov[3, 3])
rows.append(["pharmacological_sex_by_source", len(z), z.study.nunique(), beta[3], beta[3] - Z * se, beta[3] + Z * se, normal_p(beta[3] / se)])

z = d[(d.x >= 5) & d.source_group.isin(["ITP", "non-ITP pharmacological"])].copy()
female = (z.gender == "Female").astype(float).to_numpy()
non_itp = (z.source_group == "non-ITP pharmacological").astype(float).to_numpy()
x = z.x.to_numpy()
X = np.column_stack([x, x * female, x * non_itp, x * female * non_itp])
beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
se = math.sqrt(vcov[3, 3])
rows.append(["pharmacological_sex_by_source_weight_loss_ge_5pct", len(z), z.study.nunique(), beta[3], beta[3] - Z * se, beta[3] + Z * se, normal_p(beta[3] / se)])

# SD2: body-fat synthesis.
b = pd.read_excel(SOURCE, sheet_name="SD2_body_composition").copy()
x = b.bodyfat_change_points.astype(float)
y = b.lifespan_change_percent.astype(float)
study = b.pubmed_id.astype(str)
pearson_r, pearson_p = stats.pearsonr(x, y)
rows.append(["bodyfat_Pearson_r", len(b), study.nunique(), pearson_r,
             *fisher_ci(pearson_r, len(b)), pearson_p])
spearman_r, spearman_p = stats.spearmanr(x, y)
rows.append(["bodyfat_Spearman_rho", len(b), study.nunique(), spearman_r,
             np.nan, np.nan, spearman_p])
g = b.assign(study=study).groupby("study").agg(x=("bodyfat_change_points", "mean"), y=("lifespan_change_percent", "mean"))
study_r, study_p = stats.pearsonr(g.x, g.y)
rows.append(["bodyfat_study_level_r", len(b), len(g), study_r,
             *fisher_ci(study_r, len(g)), study_p])
bodyfat_study = b.pubmed_id.astype(str)
for omitted in g.index:
    subset = g.drop(index=omitted)
    result_r, result_p = stats.pearsonr(subset.x, subset.y)
    retained_arms = int((bodyfat_study != omitted).sum())
    rows.append([f"bodyfat_leave_one_study_out_{omitted}", retained_arms, len(subset),
                 result_r, *fisher_ci(result_r, len(subset)), result_p])
X = np.column_stack([np.ones(len(b)), x.to_numpy()])
beta, vcov = clustered_ols(X, y.to_numpy(), study.to_numpy())
se = math.sqrt(vcov[1, 1])
rows.append(["bodyfat_study_clustered_slope", len(b), study.nunique(), beta[1], beta[1] - Z * se, beta[1] + Z * se, normal_p(beta[1] / se)])

out = pd.DataFrame(rows, columns=["analysis", "arms", "studies", "estimate", "ci_low", "ci_high", "p_normal"])
out.to_csv(OUTPUT, index=False)
print(out.to_string(index=False))
