# CoCoNuT
 Covariate-assisted composite null hypothesis testing providing FDR control

## Installation

```R
## Install dependency packages if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")

## Install CoCoNuT
install.packages("devtools")
devtools::install_github("YanLi15/CoCoNuT")

## Load CoCoNuT
library(CoCoNuT)
```

## An numeric example

We illustrate the usage of CoCoNuT for replicability analysis of two simulated datasets incorporating auxiliary covariates.

```R
## obtained two sequence of primary p-values and two sequences of auxiliary covariates: pa, pb, x1, x2
data.obj <- rep_data_gen_multi(m = 10000, xi = c(0.9, 0.025, 0.025, 0.05), info.str =                                    c("moderate", "moderate"), mu1 = 2, mu2 = 2, mu3 = 2, 
                               mu4 = 2)
pvals1 = data.obj$pvals1
pvals2 = data.obj$pvals2
x1 = data.obj$x1
x2 = data.obj$x2

## replicability analysis exploiting the auxiliary covariates
x = cauchy(rbind(x1, x2))
res.coconut <- coconut(pvals1, pvals2, x)
coco.radj <- res.coconut$radj

## calculate the false discovery proportion and true positive proportion
states1 = data.obj$theta1
states2 = data.obj$theta2
fdr = sum(coco.radj <= 0.05 & !(states1 * states2))/ max(sum(coco.radj <= 0.05), 1)
power = sum(coco.radj <= 0.05 & (states1 * states2)) / sum(states1 * states2)
```

