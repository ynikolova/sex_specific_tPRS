# Sex-specific tPRS associations with brain structure and depressive symptoms

This repository contains R scripts examining associations of sex-specific transcriptome-based polygenic risk scores (tPRS) with structural brain measures and depressive symptoms in the Adolescent Brain Cognitive Development (ABCD) study.

The analyses include cross-sectional associations at baseline and two-year follow-up, longitudinal associations and interactions with time, non-winsorized sensitivity analyses, and selected brain structure–symptom associations and mediation models.

## Model definitions

| Model | Score adjustment | Scripts |
| --- | --- | --- |
| Model 1 | Target sex-specific tPRS, with the covariates specified in each analysis | 43 |
| Model 2 | Adds the opposite-sex tPRS | 12 |
| Model 3 | Adds the opposite-sex tPRS, original non-sex-specific tPRS (`TPRS_v1_unweighted`), and GWAS-based score (`PGC_PRS`) | 12 |

The original study README describes Models 2 and 3 as follow-up analyses of associations significant in Model 1 (FDR-adjusted p < 0.05). Some supplied scripts loop over broader sets of outcomes; the repository does not automatically restrict them to significant Model 1 findings. Consult each script's outcome selection and formulas. Additional score adjustment can occur during residualization rather than in the final mixed-model formula.

## Repository structure

| Location | Description |
| --- | --- |
| `scripts/model1/` | Primary analyses, including non-winsorized sensitivity variants |
| `scripts/model2/` | Analyses adjusted for the opposite-sex tPRS |
| `scripts/model3/` | Analyses with additional tPRS and GWAS score adjustment |
| `R/` | Shared dependency installation and path/output setup |
| `input_data/README.md` | Required private inputs and format guidance |
| `docs/script_inventory.csv` | Per-script inputs, formulas, and review flags |
| `docs/review_notes.md` | Reproducibility considerations |
| `LICENSE` | MIT license draft for code and documentation |

Each model has a `cbcl/` folder and a `freesurfer/winsorized/` folder. Model 1 also has `freesurfer/non_winsorized/`. All 67 scripts use descriptive filenames without historical date suffixes. Where names would otherwise collide, a variant suffix distinguishes them; this does not identify the final manuscript version. Run scripts independently in fresh R sessions.

## Requirements

Install R and the packages listed in `dependencies.txt`: data.table, dplyr, lme4, nlme, zoo, DescTools, and mediation. From the repository root, run:

```sh
Rscript R/install_dependencies.R
```

Exact package versions are not recorded. Refer to the review notes when reconstructing the original analysis environment.

## Data availability and preparation

Participant-level data are not included. Users must obtain the required study data through the applicable authorized access process and prepare the study-specific predictor and transformed CBCL tables.

Place the required files in `input_data/`, following the filenames and guidance in [input_data/README.md](input_data/README.md). Scripts rely on the original input column order and, in some cases, different predictor-table versions. Numbered data columns may vary with input versions, ordering, and merges, and will likely need adjustment. Verify them against `names()` on the exact table before running. Inline comments explain numeric column selections; selections whose exact names cannot be verified without the original table headers are marked accordingly.

The repository does not include the upstream score-generation or input-preparation workflow. Raw study downloads alone are therefore insufficient to run all analyses.

## Running an analysis

Open a terminal in the repository root and run a script in a fresh R session. For example:

```sh
Rscript scripts/model1/cbcl/test_tPRS_F_x_CBCL_time_1_main_by_sex.R
```

Alternatively, set the R working directory to the repository root and use:

```r
source("scripts/model1/cbcl/test_tPRS_F_x_CBCL_time_1_main_by_sex.R")
```

Results are saved as `results.csv` in `outputs/<model>/<script-name>/`. Separate output directories distinguish analysis variants. Input data, outputs, and R session files are excluded from Git tracking through `.gitignore`.

## Analysis conventions

- `MetaXcan_tPRS_F` and `MetaXcan_tPRS_M` identify score signatures. Sample sex is specified separately in each script; the supplied coding uses 1 for females and 0 for males.
- Models commonly account for age, ancestry principal components, and study site. MRI area and volume models commonly include estimated total intracranial volume (eTIV), whereas thickness models omit eTIV.
- Several longitudinal analyses use covariate-residualized outcomes. Time coding differs between CBCL and MRI analyses; consult each script.
- Scripts apply false discovery rate (FDR) correction within their specified analysis families.

## Reproducibility status

Model formulas, filtering, numeric column selections, winsorization, rounding, and FDR calculations have been retained from the analysis scripts. Comments and shared path/output setup support reuse across computers.

Static checks confirmed that the cleanup preserved executable analysis statements apart from path setup, generic file naming, and output routing. Shared setup explicitly binds `data.table::rbindlist` and `DescTools::Winsorize`, which some original scripts use without loading their packages.

## Citation

A manuscript citation and repository DOI are not yet specified. When available, we will add them here so users can cite the associated study and the version of the code used.

## License

The repository code and documentation are provided under the [MIT License](LICENSE). This permits use, modification, and redistribution, including commercial use, provided the copyright and permission notices are retained. The software is provided without warranty.

This license does not grant access to participant-level data or replace the terms governing those data. External R packages remain subject to their respective licenses.

The copyright-holder field in `LICENSE` must be completed by the repository maintainers before publication.
