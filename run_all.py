"""Run the manuscript analyses, regenerate all figures, and verify key results."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts" / "analysis"


def find_rscript(given: str | None) -> str:
    """Locate Rscript on PATH or in standard Windows installations."""
    if given:
        return given
    found = shutil.which("Rscript")
    if found:
        return found
    program_files = Path(r"C:\Program Files\R")
    candidates = sorted(program_files.glob("R-*\\bin\\Rscript.exe"), reverse=True)
    if candidates:
        return str(candidates[0])
    raise SystemExit("Rscript not found; pass --rscript with its full path")


def run(executable: str, *arguments: str | Path) -> None:
    command = [executable, *(str(argument) for argument in arguments)]
    print("\n$", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--rscript", default=None)
    args = parser.parse_args()
    python = args.python
    rscript = find_rscript(args.rscript)

    run(python, SCRIPTS / "derive_arm_level_tables.py")
    run(rscript, SCRIPTS / "analyze_control_weight_gradient.R")
    run(python, SCRIPTS / "prepare_allarm_withinarm.py")
    run(rscript, SCRIPTS / "fit_allarm_withinarm.R")
    run(rscript, SCRIPTS / "reproduce_weight_lowering_split.R")
    run(python, SCRIPTS / "prepare_predosing_quartile_analysis.py")
    run(python, SCRIPTS / "prepare_predosing_nonlowering_quartiles.py")
    run(rscript, SCRIPTS / "analyze_matched_quartiles.R")
    run(rscript, SCRIPTS / "analyze_descriptive_sensitivities.R")
    run(rscript, SCRIPTS / "analyze_compound_all13.R")
    run(rscript, SCRIPTS / "analyze_compound_common_interaction.R")
    run(python, SCRIPTS / "prepare_successful_landmarks.py")
    run(rscript, SCRIPTS / "fit_successful_landmarks.R")
    run(rscript, SCRIPTS / "analyze_comparator_cells.R")
    run(python, SCRIPTS / "reproduce_comparator_state.py")
    run(rscript, SCRIPTS / "analyze_comparator_gradient.R")
    run(rscript, SCRIPTS / "analyze_puberty_weight_moderation.R")
    run(rscript, SCRIPTS / "audit_puberty_joint.R")
    run(rscript, SCRIPTS / "fit_early_weight_change_full_cohort.R")
    run(rscript, SCRIPTS / "reproduce_maternal_adjustment.R")
    run(python, SCRIPTS / "analyze_stratum_consistency.py")
    run(python, SCRIPTS / "summarize_founder_strain_comparison.py")
    run(python, SCRIPTS / "reproduce_supplementary_analyses.py")
    run(
        rscript,
        SCRIPTS / "prepare_developmental_figure.R",
        ROOT / "data" / "inputs" / "itp_gn_earlylife_controls.csv",
        ROOT / "data" / "outputs" / "developmental_figure_estimates.csv",
    )
    run(python, SCRIPTS / "make_figures.py", "--package", ROOT,
        "--output", ROOT / "figures")
    run(python, SCRIPTS / "verify_reported_results.py")
    print("\nAll analyses and figures reproduced successfully.")


if __name__ == "__main__":
    main()
