# MoBa_data_cleaning

Files `MoBa_cleaning_MaternalHeightWeight.R` and `MoBa_imputing_MaternalHeightWeight.R` can be used as-is, on .csv files freshly made from the databases. The outputs can then be used as phenotypes or covariates themselves, or used together with `SGA-LGA_definitions_in_MoBa/` repository scripts to produce customized growth curves.

I do not take responsibility for the content/documentation of all other scripts in this repo.

Julius, 2017-05-09

## List of checks in the cleaning stage
* hard cut <130 or >210 cm height
* correct loss of meter (50 -> 150)
* correct loss of last unit (16 -> 160)
* hard cut <30 or >200 kg weight
* correct loss of last digit (9.5 -> 95)
* correct loss of decimal point (600 -> 60.0)
* remove weights with >30 kg change pre- to post-pregnancy
* remove BMI outliers on pre- and post-pregnancy weight
* remove values with large change between same mother's pregnancies

## List of steps in the imputing stage
* impute pre-/post-pregnancy weight based on the other one and KJ intake
* impute BMI based on age
* impute height as mean height
* impute weight from height and BMI
