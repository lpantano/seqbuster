pkgname <- "isomiRs"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('isomiRs')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("filter.table")
### * filter.table

flush(stderr()); flush(stdout())

### Name: filter.table
### Title: Filter isomiRs according their relative expression
### Aliases: filter.table
### Keywords: package

### ** Examples


library(isomiRs)




cleanEx()
nameEx("isomiRs-package")
### * isomiRs-package

flush(stderr()); flush(stdout())

### Name: isomiRs-package
### Title: What the package does (short line) ~~ IsomiRs analysis ~~
### Aliases: isomiRs-package isomiRs
### Keywords: package

### ** Examples


library(isomiRs)




### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
