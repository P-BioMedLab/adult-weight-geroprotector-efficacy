"""Fail fast if key regenerated outputs drift from manuscript estimands."""
from pathlib import Path
import math
import pandas as pd
root = Path(__file__).resolve().parents[2]
out = root / "data" / "outputs"
inputs = root / "data" / "inputs"
derived = {
    "weight_effects": pd.read_csv(out / "weight_effect_by_arm.csv"),
    "interactions": pd.read_csv(out / "withinarm_interactions.csv"),
    "comparator_cells": pd.read_csv(out / "comparator_effect_cells.csv"),
}
assert len(derived["weight_effects"]) == 347
assert len(derived["interactions"]) == 174
assert len(derived["comparator_cells"]) == 426
assert len(derived["comparator_cells"].query("cohort == 'C2021' and site == 'TJL'")) == 12
primary = pd.read_csv(out / "itp_allarm_primary_reproduced.csv").set_index("sex")
assert int(primary.loc["m", "arms"]) == 89 and int(primary.loc["m", "unique_mice"]) == 14301
assert int(primary.loc["f", "arms"]) == 85 and int(primary.loc["f", "unique_mice"]) == 12416
assert math.isclose(primary.loc["m", "HR"], 0.9504483787, rel_tol=0, abs_tol=1e-8)
split = pd.read_csv(out / "weight_lowering_split_reproduced.csv").set_index("analysis")
assert int(split.loc["male_weight_lowering", "arms"]) == 39
assert math.isclose(split.loc["male_weight_lowering", "HR"], 0.897457, abs_tol=1e-6)
assert int(split.loc["male_predosing_weight_lowering", "arms"]) == 27
nonlowering = split.loc["male_predosing_nonlowering"]
assert int(nonlowering.arms) == 29 and round(nonlowering.HR, 2) == 1.02
assert round(nonlowering.lo, 2) == 0.96 and round(nonlowering.hi, 2) == 1.08
assert round(nonlowering.p, 2) == 0.49
q = pd.read_csv(out / "matched_quartile_standardized_survival.csv")
assert q["gain_pct"].round(1).tolist() == [6.6, 10.1, 12.6, 20.4]
trend = pd.read_csv(out / "matched_quartile_trend.csv").iloc[0]
assert int(trend.arms) == 27 and int(trend.unique_mice) == 6080
assert math.isclose(trend.HR_per_quartile, 0.9164672687, rel_tol=0, abs_tol=1e-8)
nonlower = pd.read_csv(out / "predosing_weight_nonlowering_29arm_unique.csv")
treated_nonlower = nonlower[nonlower["Tx"].eq(1)]
assert len(treated_nonlower[["cohort", "group"]].drop_duplicates()) == 29
assert nonlower["mouse"].nunique() == 6658
comp = pd.read_csv(out / "compound_all13_descriptive.csv")
assert len(comp) == 14 and int(comp.iloc[-1].arms) == 24 and int(comp.iloc[-1].unique_mice) == 5523
assert round(comp.iloc[-1].lowest_gain_pct,1) == 5.2 and round(comp.iloc[-1].highest_gain_pct,1) == 20.9
assert int((comp.iloc[:-1].highest_gain_pct > comp.iloc[:-1].lowest_gain_pct).sum()) == 12
named = comp.set_index("compound")
for compound, lowest, highest in [
    ("Acarbose", 7.3, 28.4), ("17alpha-estradiol", 5.9, 22.9),
    ("Canagliflozin", 6.1, 23.2), ("Rapamycin", 6.7, 11.6),
]:
    assert round(named.loc[compound, "lowest_gain_pct"], 1) == lowest
    assert round(named.loc[compound, "highest_gain_pct"], 1) == highest
common = pd.read_csv(out / "compound_common_interaction.csv").iloc[0]
assert int(common.compounds) == 13 and math.isclose(common.estimate, 0.9058527, abs_tol=1e-7)
comparator = pd.read_csv(out / "comparator_state_reproduced.csv")
tertiles = comparator[comparator.split.eq(3)].set_index(["sex", "bin"])
assert round(tertiles.loc[("m", "Q1"), "mean_cell_km_median_gain_pct"], 1) == 4.1
assert round(tertiles.loc[("m", "Q3"), "mean_cell_km_median_gain_pct"], 1) == 9.3
assert round(tertiles.loc[("f", "Q1"), "mean_cell_km_median_gain_pct"], 1) == 1.6
assert round(tertiles.loc[("f", "Q3"), "mean_cell_km_median_gain_pct"], 1) == 4.1
gradient = pd.read_csv(out / "comparator_gradient_sensitivity.csv").set_index(["sex", "analysis"])
assert round(gradient.loc[("m", "within_cohort"), "p"], 3) == 0.008
assert round(gradient.loc[("m", "shared_control_clustered"), "p"], 3) == 0.078
assert round(gradient.loc[("f", "shared_control_clustered"), "p"], 2) == 0.37

age = pd.read_csv(out / "control_weight_gradient_by_age.csv").set_index(["weight_age", "sex"])
expected_age = {
    (6, "m"): (3799, 1.3919128, 1.3448974, 1.4405717, 2.3281595e-79),
    (6, "f"): (3543, 1.0771382, 1.0415789, 1.1139114, 1.4345925e-5),
    (12, "m"): (3469, 1.3575554, 1.3082173, 1.4087542, 6.5161439e-59),
    (12, "f"): (3451, 1.0915548, 1.0548034, 1.1295867, 5.3480224e-7),
    (18, "m"): (3019, 1.1221816, 1.0782903, 1.1678594, 1.4882788e-8),
    (18, "f"): (3405, 1.0495681, 1.0140265, 1.0863555, 5.9144396e-3),
    (24, "m"): (1927, 0.8437196, 0.8014716, 0.8881947, 8.9496495e-11),
    (24, "f"): (2591, 0.9764289, 0.9379186, 1.0165204, 2.4528434e-1),
}
for key, (expected_n, expected_hr, expected_lo, expected_hi, expected_p) in expected_age.items():
    assert int(age.loc[key, "n"]) == expected_n
    for column, expected in (
        ("HR_per_SD", expected_hr), ("lo", expected_lo),
        ("hi", expected_hi), ("p", expected_p),
    ):
        assert math.isclose(age.loc[key, column], expected,
                            rel_tol=1e-6, abs_tol=1e-12)
puberty = pd.read_csv(out / "itp_gn_puberty_weight_moderation.csv").set_index("sex")
assert int(puberty.loc["m", "arms"]) == 4 and int(puberty.loc["m", "unique_mice"]) == 730
assert round(puberty.loc["m", "HR"], 2) == 1.00 and round(puberty.loc["m", "p"], 2) == 0.97
assert int(puberty.loc["f", "arms"]) == 2 and int(puberty.loc["f", "unique_mice"]) == 553
assert round(puberty.loc["f", "HR"], 2) == 1.03 and round(puberty.loc["f", "p"], 2) == 0.77
joint = pd.read_csv(out / "puberty_joint_reproduced.csv").set_index("term")
assert int(joint.loc["w6z", "n"]) == 539 and round(joint.loc["w6z", "HR"], 2) == 1.46
gain = pd.read_csv(out / "early_weight_change_reproduced.csv")
male_gain = gain[(gain.sex == "m") & (gain.timing == "primary") & gain.adjusted_for_w6].iloc[0]
female_gain = gain[(gain.sex == "f") & (gain.timing == "primary") & gain.adjusted_for_w6].iloc[0]
assert int(male_gain.n_mice) == 11300 and round(male_gain.hr, 2) == 1.09
assert int(female_gain.n_mice) == 10521 and round(female_gain.hr, 2) == 1.03
assert int(male_gain.cells) == 207 and round(male_gain.ci_low, 2) == 1.06 and round(male_gain.ci_high, 2) == 1.12
assert int(female_gain.cells) == 201 and round(female_gain.ci_low, 2) == 1.01 and round(female_gain.ci_high, 2) == 1.05
maternal = pd.read_csv(out / "itp_gn_maternal_adjustment.csv")
assert (maternal.HR_weight - maternal.HR_adjusted_maternal).abs().max() < 0.01
cross = pd.read_csv(out / "supplementary_analyses_reproduced.csv").set_index("analysis")
assert int(cross.loc["pharmacological_sex_by_source", "arms"]) == 107
assert round(cross.loc["bodyfat_Pearson_r", "estimate"], 2) == -0.10
assert round(cross.loc["bodyfat_Pearson_r", "ci_low"], 2) == -0.41
assert round(cross.loc["bodyfat_Pearson_r", "ci_high"], 2) == 0.24
assert round(cross.loc["bodyfat_Pearson_r", "p_normal"], 2) == 0.58
assert round(cross.loc["bodyfat_Spearman_rho", "estimate"], 2) == 0.08
assert round(cross.loc["bodyfat_Spearman_rho", "p_normal"], 2) == 0.66
assert round(cross.loc["bodyfat_study_level_r", "estimate"], 2) == -0.21
assert round(cross.loc["bodyfat_study_level_r", "p_normal"], 2) == 0.47
expected_slopes = {
    "slope_ITP_Female": (25, 12, 0.46, 0.17, 0.75),
    "slope_ITP_Male": (31, 13, 1.03, 0.75, 1.32),
    "slope_non-ITP pharmacological_Female": (22, 15, 0.93, 0.55, 1.31),
    "slope_non-ITP pharmacological_Male": (29, 23, 0.79, 0.40, 1.18),
    "slope_dietary restriction_Female": (15, 7, 0.79, 0.62, 0.97),
    "slope_dietary restriction_Male": (12, 8, 0.81, 0.68, 0.94),
}
for analysis, (arms, studies, estimate, ci_low, ci_high) in expected_slopes.items():
    assert int(cross.loc[analysis, "arms"]) == arms
    assert int(cross.loc[analysis, "studies"]) == studies
    assert round(cross.loc[analysis, "estimate"], 2) == estimate
    assert round(cross.loc[analysis, "ci_low"], 2) == ci_low
    assert round(cross.loc[analysis, "ci_high"], 2) == ci_high
expected_medians = {
    "median_ratio_ITP_Female": 0.40,
    "median_ratio_ITP_Male": 1.27,
    "median_ratio_dietary restriction_Female": 0.67,
    "median_ratio_dietary restriction_Male": 0.74,
    "median_ratio_non-ITP pharmacological_Female": 1.00,
    "median_ratio_non-ITP pharmacological_Male": 1.29,
}
for analysis, expected in expected_medians.items():
    assert round(cross.loc[analysis, "estimate"], 2) == expected
for analysis, difference, ci_low, ci_high, p_value in [
    ("female_minus_male_ITP", -0.57, -0.98, -0.16, 0.0063),
    ("female_minus_male_non-ITP pharmacological", 0.14, -0.31, 0.60, 0.54),
    ("female_minus_male_dietary restriction", -0.02, -0.26, 0.23, 0.90),
]:
    assert round(cross.loc[analysis, "estimate"], 2) == difference
    assert round(cross.loc[analysis, "ci_low"], 2) == ci_low
    assert round(cross.loc[analysis, "ci_high"], 2) == ci_high
    assert float(f"{cross.loc[analysis, 'p_normal']:.2g}") == p_value
assert cross.index.str.startswith("bodyfat_leave_one_study_out_").sum() == 14
loo = cross[cross.index.str.startswith("bodyfat_leave_one_study_out_")]
assert loo.arms.between(30, 35).all() and loo.studies.eq(13).all()
assert (loo.p_normal > 0.05).all()
assert round(cross.loc["bodyfat_study_clustered_slope", "estimate"], 2) == -0.17
assert round(cross.loc["bodyfat_study_clustered_slope", "p_normal"], 2) == 0.57
assert round(cross.loc["pharmacological_sex_by_source", "p_normal"], 3) == 0.021
assert round(cross.loc["pharmacological_sex_by_source_weight_loss_ge_5pct", "p_normal"], 3) == 0.018

strain = pd.read_csv(out / "strain_context_summary.csv").set_index("quantity")
assert int(strain.loc["age_matched_strains", "value"]) == 27
assert round(float(strain.loc["survey_mean_weight_g", "value"]), 1) == 30.3
for quantity, expected in {
    "panel_a_strains": 27,
    "panel_a_at_or_above_itp": 6,
    "panel_b_strains": 40,
    "panel_b_at_or_above_itp": 2,
    "panel_c_strains": 40,
    "panel_c_at_or_above_itp": 1,
}.items():
    assert int(strain.loc[quantity, "value"]) == expected
founder = pd.read_csv(out / "founder_strain_weight_comparison.csv").set_index("quantity")
for quantity, expected in {
    "four_founder_mean_weight_g": 29.4,
    "itp_male_control_6m_n": 3395,
    "itp_male_control_6m_mean_weight_g": 40.4,
    "itp_percent_above_founder_mean": 37.5,
    "lowest_itp_quarter_mean_weight_g": 34.8,
    "heaviest_founder_reported_mean_weight_g": 32.0,
}.items():
    assert round(float(founder.loc[quantity, "value"]), 1) == expected
print("Verified the principal, class, quartile, compound, comparator, developmental and cross-study results.")
