# PAthogenic VAriant Reporter (PAVAR)

PAVAR classifies pathogenic rare variants based on the guideline of ACMG/AMP and estimate potential carrier in U.S. population based on Hardy-Weinberg equilibrium (HWE).

The program is written in python except the HWE implemented in R.

- Prerequisite:
  - [InterVar](https://github.com/WGLab/InterVar)
    - Installation steps:
  - R library dependencies: the hwe.R script automatically installes the packages if not installed.
    - [HardyWeinberg](https://cran.r-project.org/web/packages/HardyWeinberg/index.html)
    - [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
    - [data.table](https://cran.r-project.org/web/packages/data.table/)

- Now you are ready to run the [paver.py](https://github.com/dauss75/pavar/blob/master/bin/paver.py) script.
