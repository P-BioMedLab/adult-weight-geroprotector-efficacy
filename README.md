# Adult body weight and geroprotector efficacy

Code and data supporting the analyses in *Adult body weight modifies geroprotector efficacy in genetically heterogeneous mice*.

## Reproduction

Tested with Python 3.13, R 4.5.1 and the R package `survival` 3.8-3.

```bash
python -m pip install -r requirements.txt
python run_all.py
```

The command derives the arm-level tables from the harmonized mouse data, regenerates the reported analyses and all four figures, and runs the numerical checks. Use `--python` or `--rscript` to specify either executable. To regenerate only the figures after running the analyses:

```bash
python scripts/analysis/make_figures.py --package . --output figures
```

## Contents

- `data/inputs`: harmonized mouse-level and external source data
- `data/outputs`: regenerated arm-level tables, estimates and tabulations
- `scripts/analysis`: analysis code; `make_figures.py` produces all four figures
- `figures`: manuscript figures in PNG and TIFF formats
- `supplementary/Supplementary_Data_1-6.xlsx`: Supplementary Data 1-6
- `ANALYSIS_HIERARCHY.md`: inferential hierarchy recorded at repository freeze
- `checksums.sha256`: SHA-256 manifest for the archived package

Primary ITP data are available from the [Mouse Phenome Database ITP portal](https://phenome.jax.org/centers/ITP) (projects ITP1 and ITP2004-ITP2021) and the [GeneNetwork HET3-ITP resource](https://genenetwork.org/show_trait?dataset=HET3-ITPPublish&trait_id=50601). Founder-strain data are from [Mouse Phenome Database project Reed1](https://phenome.jax.org/projects/Reed1); archived strain-survey inputs allow offline figure generation.

Copyright © 2026 the authors. Code is licensed under AGPL-3.0-only. Source data remain subject to the terms of their respective providers; see `DATA_LICENSE.md`.
