
# **ULV**

ULV is a two-stage robust statistical framework that incorporates U statistic with a latent variable model for the purpose of comparing the difference of clustered data. 
In the first stage, ULV utilizes a nonparametric rank-based method to calculate the between-subject difference between individuals from distinct groups. 
In the second stage, ULV postulates a latent level for each subject via a parametric model. 
The analysis will benefit from ULV because the rank-based method in the first stage is robust to non-normal distribution of data, while the latent variable model in the second stage accounts for correlations in the between-subject difference, thereby characterize the data dependence caused by clustering.

## Installation

You can install the current version of ULV from GitHub with:
```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yu-zhaoxia/ULV")
```

## Example

This is a basic example which shows you how to apply ULV to a single-cell RNA-seq data:
```r
library(ULV)
data('example_data')
count = example_data$count_matrix
meta = example_data$metadata

res_table = fit_ULV(count, meta, normalize=TRUE, 
                subject_name = 'donor', cond_name = 'group_per_sample', 
                ctrl_cond = 'mild', case_cond = 'severe', 
                weighted = TRUE, covariate_name_list=c('age_yr','sex'))
```
