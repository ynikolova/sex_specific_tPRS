# Review notes

This cleanup preserves the original statistical calculations. These observations are review flags, not confirmed result errors. The per-script inventory identifies affected scripts.

- **FDR precision:** many scripts round raw p-values to three decimals before applying FDR. Consider adjusting full-precision p-values, then rounding only for display. This would change results and was not performed.
- **Sample counts:** many reported n values subtract missingness in one column rather than count complete cases actually used by the fitted model. Some longitudinal counts use rows, others approximate participant counts. Verify against fitted-model observations and unique participants.
- **Random effects:** `nlme::lme` interprets a named list of grouping factors as a hierarchy. Review `list(SID=~time,site=~1)` against the intended site/participant hierarchy; it is not a general declaration of crossed effects.
- **Column positions:** many selections and winsorization ranges use numeric indices. Changes in predictor-table versions or column order can select different variables silently. The 2023 and 2026 predictor tables are both required by different variants.
- **Events and time:** several MRI longitudinal scripts classify every non-baseline event as follow-up and retain any duplicated SID. Confirm that only the two intended events are present. CBCL commonly uses 1/2 while MRI commonly uses 0/1.
- **Coefficient extraction:** fixed-effect rows are selected by number; confirm they correspond to the intended term. Reported t and p are often rounded before storage.
- **Mediation:** `a_t` and `b_t` contain coefficient estimates (column 1), not t statistics. The standalone `fit.b` is unadjusted for tPRS, whereas the mediation outcome model includes tPRS. Bootstrap seeds are not set. Statistical mediation alone does not establish causality.
- **Variants and naming:** all versions were retained. The MRI filename containing `separate_T1_T1` may be a typo; the analysis label was retained. Some non-winsorized scripts retain winsorized phenotype labels or inactive Winsorize lines. Do not select a manuscript version solely from the filename.
- **Dependencies:** no original sessionInfo or lockfile was supplied. DescTools Winsorize uses default quantile limits; check package behavior and missing-data handling in the original environment. Official reference: https://andrisignorell.github.io/DescTools/reference/Winsorize.html
- **Input provenance:** transformed CBCL tables and predictor preparation were not supplied. Exact source schemas and score construction cannot be reconstructed from these scripts alone.

No inferential fixes were applied during this cleanup. Full reproducibility requires the matching private input tables, original package versions, and an execution check.

## Numbered-column annotations

Inline comments now identify the purpose of each numeric column selection, including predictor/PC blocks, MRI conversion ranges, merged-table winsorization offsets, and posthoc difference columns. Exact names unavailable without private input headers are explicitly marked unverified. Use `names(table)[indices]` to inspect these locally. Matrix positions used for coefficients and confidence intervals are also explained. No computer-specific absolute paths are present in the scripts; relative repository paths remain necessary for inputs and shared setup.

## Combined models

- Model definitions derive from the supplied study README. That README contained computer-specific paths and was replaced with the shared public README; unavailable upstream calculation scripts are described without local paths.
- The Model 2 archive also contained byte-identical copies of Models 1 and 3. One copy of each script is included. Model 1 retains the previously cleaned version, and Model 3 uses its separately supplied archive.
- All 24 additional scripts retain original calculations; file references use generic names. Shared setup binds rbindlist and Winsorize explicitly because several scripts omit the package declarations needed for these calls in a fresh session.
- Some adjusted models include additional scores in residualization rather than the final model. Mediation-model formulas may therefore remain identical between models while the residualized inputs differ.
- Review whether the outcome loops implement the intended restriction to significant Model 1 findings. No significance-based gating was added.
- Across-model comparison requires care when additional covariates change complete-case samples. Statistical calculations were not revised.

## Adjust numbered columns for your inputs

Every executable numbered data-column selection now has an adjacent portability notice: indices reflect the original layout, may vary by release/version/order/merge, and will likely need adjustment. Check the intended fields in the preceding comments against `names()` on the exact table being indexed. Numeric coefficient rows also require checking after model changes; confidence-interval indices represent lower/upper bounds, not input-data columns. Description-row removal must be checked against the actual export format.

## Generic file names

Inputs use symptoms_baseline.csv, symptoms_followup.csv, predictors_primary.csv, predictors_sensitivity.csv, and structural_measures.txt. The two predictor versions remain distinct because their schemas may differ. Rename your local input files accordingly; expected contents and column assumptions are unchanged. Each script writes results.csv in its own model/analysis output directory. Script date stamps were removed; colliding names receive a variant suffix.
