# Private inputs

Place the following files here locally. No participant data are included.

* `symptoms\_followup.csv`
* `symptoms\_baseline.csv`
* `predictors\_primary.csv`
* `predictors\_sensitivity.csv`
* `structural\_measures.txt`

These are study-specific processed tables, not interchangeable raw downloads. Preserve original column order: scripts select many columns by number. Predictor preparation, score generation, and CBCL transformation code were not included. Inspect each script and the inventory for required covariates and fields.

Common fields include SID, sex (0 male, 1 female as coded in the scripts), has\_rel, site, interview\_age, eTIV, eigenvec.1 through eigenvec.10, MetaXcan\_tPRS\_F, and MetaXcan\_tPRS\_M. FreeSurfer input expects src\_subject\_id, eventname, smri\_\* measures and a descriptive first row that the code removes. Exact column order and additional variables must match the original tables.

Models 2 and 3 additionally require the opposite-sex MetaXcan tPRS in the relevant predictor/CBCL tables. Model 3 also requires TPRS\_v1\_unweighted and PGC\_PRS.

