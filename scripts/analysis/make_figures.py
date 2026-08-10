"""Regenerate all four manuscript figures from packaged analysis outputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BLUE = "#2F7ED8"
ORANGE = "#F26B38"
GREY = "#8F8C85"
INK = "#1A1A1A"
GRID = "#E3E3E3"
SLATE = "#55534E"
PRIMARY_BLUE = "#0072B2"
SENSITIVITY_GREY = "#A7A9AC"
DEVELOPMENTAL_BLUE = "#1F5FA9"
DEVELOPMENTAL_ORANGE = "#D4692A"
FOUNDERS = ["BALB/cByJ", "C57BL/6J", "C3H/HeJ", "DBA/2J"]


def save_figure(fig, output: Path, stem: str, dpi: int = 300) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{stem}.tif", dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def base_style(font_size: float = 8.7) -> None:
    plt.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.family": "DejaVu Sans", "font.size": font_size,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": True, "axes.spines.bottom": True,
        "axes.edgecolor": "#AAA79E", "axes.labelcolor": SLATE,
        "xtick.color": "#77746E", "ytick.color": "#77746E",
        "xtick.major.size": 0, "ytick.major.size": 0,
        "xtick.minor.size": 0, "ytick.minor.size": 0,
    })


def figure1(package: Path, output: Path) -> None:
    """All-arm treatment-by-weight estimates."""
    data = package / "data"
    interactions = pd.read_csv(data / "outputs" / "withinarm_interactions.csv")
    weight = pd.read_csv(data / "outputs" / "weight_effect_by_arm.csv")
    post = weight[weight.post_dosing.astype(str).str.lower().eq("true")].copy()
    post["weighted_diff"] = post.diff_g * post.n_tx
    lowering = post.groupby(["cohort", "group"], as_index=False).agg(
        weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum")
    )
    lowering["weight_lowering"] = lowering.weighted_diff / lowering.n_tx < -1
    interactions = interactions.merge(
        lowering[["cohort", "group", "weight_lowering"]],
        on=["cohort", "group"], how="left", validate="many_to_one"
    )
    interactions["weight_lowering"] = interactions.weight_lowering.eq(True)
    observed = interactions.groupby("sex").weight_lowering.sum().to_dict()
    if observed != {"m": 39, "f": 37}:
        raise AssertionError(f"Unexpected weight-lowering counts: {observed}")
    pooled = pd.read_csv(data / "outputs" / "itp_allarm_primary_reproduced.csv").set_index("sex")

    base_style(10)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 6.2), sharex=True)
    panels = zip(axes, ["m", "f"], ["a", "b"],
                 ["Males (89 arms)", "Females (85 arms)"])
    for ax, sex, letter, title in panels:
        frame = interactions[interactions.sex == sex].sort_values("int_hr").reset_index(drop=True)
        y = np.arange(len(frame))
        colours = np.where(frame.weight_lowering, BLUE, GREY)
        for yi, lo, hi, colour in zip(y, frame.int_lo, frame.int_hi, colours):
            ax.plot([max(lo, 0.35), min(hi, 2.1)], [yi, yi], color=colour, alpha=0.25, lw=0.8)
        ax.scatter(frame.int_hr.clip(0.35, 2.1), y, c=colours, s=12, alpha=0.8, edgecolors="none")
        ax.axvline(1, color="#AAA79E", lw=1.2)
        estimate = pooled.loc[sex]
        ax.errorbar(estimate.HR, -7,
                    xerr=[[estimate.HR - estimate.lo], [estimate.hi - estimate.HR]],
                    fmt="D", ms=8, color=BLUE, capsize=3, lw=2)
        ax.text(0.37, -11, f"Pooled HR {estimate.HR:.2f} ({estimate.lo:.2f}–{estimate.hi:.2f})",
                fontsize=9, fontweight="bold")
        ax.set_ylim(len(frame) + 2, -15)
        ax.set_yticks([])
        ax.set_title(letter, loc="left", fontweight="bold", fontsize=13, pad=22)
        ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 5),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=11, color="#333333")
        ax.set_xlabel("Treatment × weight interaction HR per +1 s.d.")
        ax.grid(axis="x", color="#E4E1DA", lw=0.7)
    fig.tight_layout(rect=[0, 0.04, 1, 0.98])
    save_figure(fig, output, "Figure_1_all_arm_interactions")


def figure2(package: Path, output: Path) -> None:
    """Matched pre-treatment quartiles and pooled sensitivity."""
    data = package / "data" / "outputs"
    frame = pd.read_csv(data / "matched_quartile_standardized_survival.csv")
    trend = pd.read_csv(data / "matched_quartile_trend.csv").iloc[0]
    sensitivity = pd.read_csv(data / "descriptive_sensitivities_reproduced.csv")
    conventional = (sensitivity[
        sensitivity.analysis.eq("weight-lowering conventional quartiles")
    ].sort_values("bin").gain_pct.to_numpy())
    if conventional.shape != (4,):
        raise AssertionError("Expected four conventional-quartile sensitivity estimates")
    x = np.arange(4)
    labels = ["Lightest\nquarter", "Second\nquarter", "Third\nquarter", "Heaviest\nquarter"]

    base_style(8.7)
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 4.10), dpi=300,
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    # Match the male/female palette used in Figure 4 while retaining the
    # distinct blue/grey encoding for the standardization comparison in panel b.
    ax.plot(x, frame.standardized_control_median, "o-", lw=2.8, ms=8,
            color=DEVELOPMENTAL_BLUE, label="Control")
    ax.plot(x, frame.standardized_treated_median, "o-", lw=2.8, ms=8,
            color=DEVELOPMENTAL_ORANGE, label="Treated")
    for i, row in frame.iterrows():
        first = i == 0
        ax.annotate(f"{row.mean_control_weight:.1f} g", (i, row.standardized_control_median),
                    xytext=(6 if first else 0, -19), textcoords="offset points",
                    ha="left" if first else "center", color=SLATE, fontsize=8.2)
        ax.annotate(f"+{row.gain_pct:.1f}%",
                    (i, (row.standardized_control_median + row.standardized_treated_median) / 2),
                    xytext=(7, 0), textcoords="offset points", va="center",
                    color="#6A401E", fontsize=8.5, fontweight="bold")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Standardized median lifespan (days)")
    ax.set_ylim(620, 990)
    ax.set_yticks(range(650, 951, 50))
    ax.grid(axis="y", color="#DEDBD3", lw=0.8)
    ax.tick_params(axis="both", labelsize=8.6)
    ax.legend(frameon=False, fontsize=8.2)
    ax.set_title("a", loc="left", fontsize=11, fontweight="bold", pad=12)
    ax.text(0, 1.025, "Matched quartile survival", transform=ax.transAxes,
            fontsize=8.5, color="#333333")

    ax = axes[1]
    width = 0.36
    ax.bar(x - width / 2, frame.gain_pct, width=width, color=PRIMARY_BLUE, label="Standardized")
    ax.bar(x + width / 2, conventional, width=width, color=SENSITIVITY_GREY, label="Pooled sensitivity")
    ax.set_xticks(x, ["Lowest", "Second", "Third", "Highest"])
    ax.tick_params(axis="both", labelsize=8.6)
    ax.set_ylabel("Median lifespan gain (%)")
    ax.set_ylim(0, 24)
    ax.grid(axis="y", color="#DEDBD3", lw=0.8)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2,
              fontsize=8.2, frameon=False, handlelength=1.5, columnspacing=1.3)
    ax.set_title("b", loc="left", fontsize=11, fontweight="bold", pad=12)
    ax.text(0, 1.025, "Robust to standardization", transform=ax.transAxes,
            fontsize=8.5, color="#333333")
    ax.text(0.02, 0.96,
            f"Trend HR {trend.HR_per_quartile:.2f} ({trend.lo:.2f}–{trend.hi:.2f})\nP = {trend.p:.4f}",
            transform=ax.transAxes, va="top", fontsize=8.2, fontweight="bold")
    fig.tight_layout(rect=[0, 0.12, 1, 0.93])
    save_figure(fig, output, "Figure_2_matched_quartiles")


def load_strain_means(path: Path, sex: str) -> dict[str, float]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    rows = raw.get("strainmeans", raw) if isinstance(raw, dict) else raw
    return {row["strain"]: float(row["mean"]) for row in rows
            if str(row.get("sex", "")).lower().startswith(sex) and row.get("mean") is not None}


def figure3(package: Path, output: Path) -> None:
    """Laboratory-strain body-weight and adiposity context."""
    inputs = package / "data" / "inputs"
    weight = load_strain_means(inputs / "mpd_32603_strainmeans.json", "m")
    age = load_strain_means(inputs / "mpd_32680_strainmeans.json", "m")
    matched = {strain: value for strain, value in weight.items()
               if strain in age and age[strain] <= 36.0}
    if len(matched) != 27 or not all(founder in matched for founder in FOUNDERS):
        raise AssertionError(f"Expected 27 age-matched strains including all founders, found {len(matched)}")
    survey_mean = float(np.mean(list(matched.values())))
    if round(survey_mean, 1) != 30.3:
        raise AssertionError(f"Unexpected survey mean: {survey_mean}")

    quartiles = pd.read_csv(package / "data" / "outputs" /
                            "matched_quartile_standardized_survival.csv")
    quartile_marks = list(zip(
        ["Lightest quarter", "Second", "Third", "Heaviest quarter"],
        quartiles.sort_values("quartile").mean_control_weight.astype(float),
    ))
    body = pd.read_excel(package / "supplementary" / "Supplementary_Data_1-6.xlsx",
                         sheet_name="SD2_body_composition")
    rapa = body[(body["Intervention"] == "Rapamycin") &
                (body["strain"] == "UM-HET3") &
                body["dosage"].astype(str).str.contains("42 ppm", regex=False)]
    adiposity = rapa.set_index("gender")["bodyfat_pct_control"].astype(float)
    if not {"Male", "Female"}.issubset(adiposity.index):
        raise AssertionError("Missing UM-HET3 control adiposity values")

    panels = [
        ("a", "Body weight, males", matched, "Body weight (g)",
         quartile_marks,
         survey_mean),
        ("b", "Body fat, males", load_strain_means(inputs / "mpd_10331_strainmeans.json", "m"),
         "Body fat (%)", [("ITP UM-HET3", adiposity.loc["Male"])], None),
        ("c", "Body fat, females", load_strain_means(inputs / "mpd_10331_strainmeans.json", "f"),
         "Body fat (%)", [("ITP UM-HET3", adiposity.loc["Female"])], None),
    ]

    base_style(10)
    plt.rcParams["axes.spines.left"] = False
    # Preserve the manuscript image aspect ratio and pixel dimensions at 220 dpi.
    fig, axes = plt.subplots(3, 1, figsize=(2344 / 220, 1755 / 220))
    for ax, (letter, title, values_by_strain, x_label, marks, reference_mean) in zip(axes, panels):
        values = np.array(list(values_by_strain.values()))
        names = list(values_by_strain.keys())
        y = np.random.default_rng(0).uniform(-0.17, 0.17, len(values))
        founder_mask = np.array([name in FOUNDERS for name in names])
        ax.scatter(values[~founder_mask], y[~founder_mask], s=26, c=GREY, alpha=0.55,
                   edgecolors="none", label=f"other strains ({(~founder_mask).sum()})", zorder=2)
        ax.scatter(values[founder_mask], y[founder_mask], s=64, c=BLUE,
                   edgecolors="white", linewidths=0.8, label="UM-HET3 founders (4)", zorder=4)
        for name, value, yy in sorted(zip(names, values, y), key=lambda item: item[1]):
            if name in FOUNDERS:
                below = FOUNDERS.index(name) % 2 == 0
                ax.annotate(name, (value, yy), xytext=(0, -17 if below else 13),
                            textcoords="offset points", ha="center", fontsize=7.2,
                            color="#3A5F8A", zorder=6)
        if reference_mean is not None:
            ax.plot([reference_mean, reference_mean], [-0.46, 0.24], ls=(0, (4, 3)),
                    color=SLATE, lw=1.6, zorder=5)
            ax.annotate(f"Survey mean {reference_mean:.1f} g", (reference_mean, -0.47),
                        ha="center", va="top", fontsize=7.6, color=SLATE,
                        fontweight="bold", zorder=7)
        for index, (label, value) in enumerate(marks):
            high = 0.60 if index % 2 == 0 else 0.92
            ax.axvline(value, color=ORANGE, lw=2.0, ymin=0.34, ymax=high, zorder=5)
            ax.annotate(f"{label}\n{value:.1f}", (value, high * 1.55 - 0.55), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=7.6, color="#A8471C", fontweight="bold", zorder=7)
        ax.set_ylim(-0.95 if reference_mean is not None else -0.72, 1.10)
        ax.set_yticks([])
        ax.set_xlabel(x_label)
        ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=11.5)
        ax.grid(axis="x", color="#E4E1DA", lw=0.7, zorder=0)
        if len(marks) == 1:
            note = f"{int((values >= marks[0][1]).sum())} of {len(values)} strains at or above the ITP value"
        else:
            note = (f"{int((values >= marks[0][1]).sum())} of {len(values)} strains reach the lightest "
                    f"ITP quarter, none reaches the heaviest. The lightest quarter is "
                    f"{100 * (marks[0][1] / reference_mean - 1):.0f}% above the survey mean")
        ax.annotate(note, (0.995, 0.02), xycoords="axes fraction", ha="right",
                    fontsize=8.0, color=SLATE)
        if letter == "a":
            ax.legend(frameon=False, fontsize=8, loc="upper left",
                      bbox_to_anchor=(0, 1.02), handletextpad=0.3)
    fig.tight_layout(rect=[0, 0.01, 1, 0.99])
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "Figure_3_strain_context.png", dpi=220, facecolor="white")
    fig.savefig(output / "Figure_3_strain_context.tif", dpi=600, facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def figure4(package: Path, output: Path) -> None:
    """Developmental decomposition in order of first citation."""
    data = package / "data" / "outputs"
    age = pd.read_csv(data / "developmental_figure_estimates.csv")
    control = pd.read_csv(data / "control_weight_gradient_by_age.csv")
    moderation_early = pd.read_csv(data / "itp_gn_puberty_weight_moderation.csv")
    moderation_adult = pd.read_csv(data / "successful_landmarks_reproduced.csv")

    early_age = age[age.panel.eq("age")].set_index("sex")
    control_6_12 = control[control.weight_age.isin([6, 12])]
    gradients = {
        ("m", 1.4): early_age.loc["m", ["HR", "lo", "hi"]].astype(float).to_numpy(),
        ("f", 1.4): early_age.loc["f", ["HR", "lo", "hi"]].astype(float).to_numpy(),
    }
    for row in control_6_12.itertuples(index=False):
        gradients[(row.sex, float(row.weight_age))] = np.array([row.HR_per_SD, row.lo, row.hi], float)

    decomp = age[age.panel.eq("decomposition")].copy()
    expected_labels = ["1.4-month alone", "6-month alone", "1.4-month adjusted", "6-month adjusted"]
    if decomp.label.tolist() != expected_labels:
        raise AssertionError(f"Unexpected decomposition rows: {decomp.label.tolist()}")

    adult12 = moderation_adult[(moderation_adult.analysis == "published extenders") &
                               (moderation_adult.landmark == 12)].set_index("sex")
    early_mod = moderation_early.set_index("sex")
    moderation_rows = [
        ("1.4-mo\nmale", *early_mod.loc["m", ["HR", "lo", "hi"]].astype(float), DEVELOPMENTAL_BLUE),
        ("12-mo\nmale", *adult12.loc["m", ["HR", "lo", "hi"]].astype(float), DEVELOPMENTAL_BLUE),
        ("1.4-mo\nfemale", *early_mod.loc["f", ["HR", "lo", "hi"]].astype(float), DEVELOPMENTAL_ORANGE),
        ("12-mo\nfemale", *adult12.loc["f", ["HR", "lo", "hi"]].astype(float), DEVELOPMENTAL_ORANGE),
    ]

    base_style(8.5)
    fig = plt.figure(figsize=(7.09, 2.95))
    grid = fig.add_gridspec(1, 3, wspace=0.44, left=0.075, right=0.985, top=0.85, bottom=0.235)

    ax = fig.add_subplot(grid[0, 0])
    ages = [1.4, 6, 12]
    male = np.array([gradients[("m", value)] for value in ages])
    female = np.array([gradients[("f", value)] for value in ages])
    ax.fill_between(ages, male[:, 1], male[:, 2], color=DEVELOPMENTAL_BLUE, alpha=0.15, lw=0)
    ax.fill_between(ages, female[:, 1], female[:, 2], color=DEVELOPMENTAL_ORANGE, alpha=0.15, lw=0)
    ax.plot(ages, male[:, 0], "-o", color=DEVELOPMENTAL_BLUE, lw=1.7, ms=4.4,
            mfc=DEVELOPMENTAL_BLUE, mec="white", mew=1)
    ax.plot(ages, female[:, 0], "-^", color=DEVELOPMENTAL_ORANGE, lw=1.7, ms=4.4,
            mfc=DEVELOPMENTAL_ORANGE, mec="white", mew=1)
    ax.axhline(1, color="#999999", lw=0.9)
    ax.set_xticks(ages, ["1.4", "6", "12"])
    ax.set_xlim(0.2, 13.2)
    ax.set_ylim(0.86, 1.50)
    ax.set_xlabel("Age at weighing (months)", fontsize=8.8, color=INK)
    ax.set_ylabel("Mortality hazard per +1 SD", fontsize=8.8, color=INK)
    ax.annotate("males", (12, male[-1, 0]), xytext=(-4, 10), textcoords="offset points",
                color=DEVELOPMENTAL_BLUE, fontsize=7.8, fontweight="bold", ha="right")
    ax.annotate("females", (12, female[-1, 0]), xytext=(-4, -15), textcoords="offset points",
                color=DEVELOPMENTAL_ORANGE, fontsize=7.8, fontweight="bold", ha="right")
    ax.set_title("a", loc="left", fontsize=10.5, fontweight="bold", pad=11)
    ax.text(0, 1.035, "Weight predicts mortality", transform=ax.transAxes, fontsize=8.2, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.6)

    ax = fig.add_subplot(grid[0, 1])
    y = np.arange(len(moderation_rows))[::-1]
    for yi, (_, hr, lo, hi, colour) in zip(y, moderation_rows):
        ax.plot([lo, hi], [yi, yi], color=colour, lw=1.6)
        ax.plot(hr, yi, "o", color=colour, ms=6.4, mec="white", mew=1.2, zorder=3)
    ax.axvline(1, color="#999999", lw=0.9)
    ax.set_yticks(y, [row[0] for row in moderation_rows], fontsize=6.8)
    ax.set_ylim(-0.6, len(moderation_rows) - 0.4)
    ax.set_xlim(0.80, 1.28)
    ax.set_xticks([0.9, 1.0, 1.1, 1.2])
    ax.set_xlabel("Treatment × weight HR", fontsize=8.8, color=INK)
    ax.set_title("b", loc="left", fontsize=10.5, fontweight="bold", pad=11)
    ax.text(0, 1.035, "Only adult weight modifies benefit", transform=ax.transAxes,
            fontsize=8.2, color=INK)
    ax.grid(axis="x", color=GRID, lw=0.6)

    ax = fig.add_subplot(grid[0, 2])
    labels = ["1.4-mo\nalone", "6-mo\nalone", "1.4-mo\nadj.", "6-mo\nadj."]
    colours = [GREY, DEVELOPMENTAL_BLUE, GREY, DEVELOPMENTAL_BLUE]
    x = np.arange(4)
    for index, row in enumerate(decomp.itertuples(index=False)):
        ax.plot([x[index], x[index]], [row.lo, row.hi], color=colours[index], lw=1.6)
        ax.plot(x[index], row.HR, "o", color=colours[index], ms=6.4,
                mec="white", mew=1.2, zorder=3)
    ax.axhline(1, color="#999999", lw=0.9)
    ax.axvline(1.5, color="#CCCCCC", lw=0.9, ls=":")
    ax.set_xticks(x, labels, fontsize=6.8)
    ax.set_ylabel("Mortality hazard per +1 SD", fontsize=8.8, color=INK)
    ax.set_ylim(0.82, 1.78)
    ax.text(0.5, 1.70, "separate", ha="center", fontsize=7.6, color=INK)
    ax.text(2.5, 1.70, "adjusted", ha="center", fontsize=7.6, color=INK)
    ax.set_title("c", loc="left", fontsize=10.5, fontweight="bold", pad=11)
    ax.text(0, 1.035, "Adult, not early-life, weight", transform=ax.transAxes,
            fontsize=8.2, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.6)

    for axis in fig.axes:
        # Use one tick-label scale across all three panels. Previously panel a
        # inherited 8.5 pt ticks while panel b's y labels and panel c's x
        # labels were 6.8 pt, which made the composite look mismatched.
        axis.tick_params(axis="both", which="major", labelsize=7.8,
                         colors=INK, labelcolor=INK)
        axis.xaxis.label.set_color(INK)
        axis.yaxis.label.set_color(INK)
        axis.set_axisbelow(True)
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    save_figure(fig, output, "Figure_4_developmental", dpi=600)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    package = args.package.resolve()
    output = (args.output or package / "figures").resolve()
    figure1(package, output)
    figure2(package, output)
    figure3(package, output)
    figure4(package, output)
    print(f"Regenerated four manuscript figures in {output}")


if __name__ == "__main__":
    main()
