pkgname <- "emdi"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('emdi')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("benchmark")
### * benchmark

flush(stderr()); flush(stdout())

### Name: benchmark
### Title: Benchmark function
### Aliases: benchmark

### ** Examples

# Loading data - population and sample data
data("eusilcA_popAgg")
data("eusilcA_smpAgg")

# Combine sample and population data
combined_data <- combine_data(pop_data = eusilcA_popAgg, pop_domains = "Domain",
                             smp_data = eusilcA_smpAgg, smp_domains = "Domain")

# Estimate Fay-Herriot model
fh_std <- fh(fixed = Mean ~ cash + self_empl, vardir = "Var_Mean",
combined_data = combined_data, domains = "Domain", method = "ml", 
MSE = TRUE)  

# Benchmark the point estimates

# Example 1: Receive data frame with point estimates and their benchmarked results
fh_bench <- benchmark(fh_std, benchmark = 20140.09, 
share = eusilcA_popAgg$ratio_n, type = "ratio")

# Example 2: Add benchmarked results to fh object
fh_bench <- benchmark(fh_std, benchmark = 20140.09, 
share = eusilcA_popAgg$ratio_n, type = "ratio", overwrite = TRUE)



cleanEx()
nameEx("compare_plot.emdi")
### * compare_plot.emdi

flush(stderr()); flush(stdout())

### Name: compare_plot.emdi
### Title: Shows plots for the comparison of estimates
### Aliases: compare_plot.emdi

### ** Examples




cleanEx()
nameEx("data_transformation")
### * data_transformation

flush(stderr()); flush(stdout())

### Name: data_transformation
### Title: Tranforms dependent variables
### Aliases: data_transformation

### ** Examples

# Loading data - sample data
data("eusilcA_smp")

# Transform dependent variable in sample data with Box-Cox transformation
transform_data <- data_transformation(eqIncome ~ gender + eqsize + cash + 
self_empl + unempl_ben + age_ben + surv_ben + sick_ben + dis_ben + rent + 
fam_allow + house_allow + cap_inv + tax_adj, eusilcA_smp, "box.cox", 0.7)



cleanEx()
nameEx("direct")
### * direct

flush(stderr()); flush(stdout())

### Name: direct
### Title: Direct estimation of disaggregated indicators
### Aliases: direct

### ** Examples




cleanEx()
nameEx("ebp")
### * ebp

flush(stderr()); flush(stdout())

### Name: ebp
### Title: Empirical Best Prediction for disaggregated indicators
### Aliases: ebp

### ** Examples




cleanEx()
nameEx("estimators.emdi")
### * estimators.emdi

flush(stderr()); flush(stdout())

### Name: estimators.emdi
### Title: Presents point, MSE and/or CV estimates of an emdiObject
### Aliases: estimators.emdi

### ** Examples




cleanEx()
nameEx("fh")
### * fh

flush(stderr()); flush(stdout())

### Name: fh
### Title: Standard and extended Fay-Herriot models for disaggregated
###   indicators
### Aliases: fh

### ** Examples




cleanEx()
nameEx("head.estimators.emdi")
### * head.estimators.emdi

flush(stderr()); flush(stdout())

### Name: head.estimators.emdi
### Title: Returns the first part of predicted indicators and, if chosen,
###   of MSE and CV estimators.
### Aliases: head.estimators.emdi

### ** Examples




cleanEx()
nameEx("map_plot")
### * map_plot

flush(stderr()); flush(stdout())

### Name: map_plot
### Title: Visualizes regional disaggregated estimates on a map
### Aliases: map_plot

### ** Examples




cleanEx()
nameEx("plot.emdi")
### * plot.emdi

flush(stderr()); flush(stdout())

### Name: plot.emdi
### Title: Plots for an emdi object
### Aliases: plot.emdi

### ** Examples




cleanEx()
nameEx("spatialcor.tests")
### * spatialcor.tests

flush(stderr()); flush(stdout())

### Name: spatialcor.tests
### Title: Spatial autocorrelation tests
### Aliases: spatialcor.tests

### ** Examples

# Loading data - sample data and proximity matrix
data("eusilcA_smpAgg")
data("eusilcA_prox")

# Compute spatial correlation tests 
spatialcor.tests(direct = eusilcA_smpAgg$Mean, 
corMatrix = eusilcA_prox)



cleanEx()
nameEx("step.fh")
### * step.fh

flush(stderr()); flush(stdout())

### Name: step.fh
### Title: Method 'step.fh' selects a Fay-Herriot model by different
###   information criteria in a stepwise algorithm.
### Aliases: step.fh

### ** Examples




cleanEx()
nameEx("subset.estimators.emdi")
### * subset.estimators.emdi

flush(stderr()); flush(stdout())

### Name: subset.estimators.emdi
### Title: Subsets an estimators.emdi object
### Aliases: subset.estimators.emdi

### ** Examples




cleanEx()
nameEx("summary.emdi")
### * summary.emdi

flush(stderr()); flush(stdout())

### Name: summary.emdi
### Title: Summarizes an emdiObject
### Aliases: summary.emdi

### ** Examples




cleanEx()
nameEx("tail.estimators.emdi")
### * tail.estimators.emdi

flush(stderr()); flush(stdout())

### Name: tail.estimators.emdi
### Title: Returns the last part of predicted indicators and, if chosen, of
###   MSE and CV estimators.
### Aliases: tail.estimators.emdi

### ** Examples




cleanEx()
nameEx("write.excel")
### * write.excel

flush(stderr()); flush(stdout())

### Name: write.excel
### Title: Exports an emdiObject to an Excel file or OpenDocument
###   Spreadsheet
### Aliases: write.excel write.ods

### ** Examples





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
