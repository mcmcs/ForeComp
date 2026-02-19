pkgname <- "ForeComp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ForeComp-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ForeComp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Plot_Tradeoff")
### * Plot_Tradeoff

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Plot_Tradeoff
### Title: Visualizes the size distortion maximum power loss tradeoff from
###   the Diebold-Mariano test for equal predictive accuracy
### Aliases: Plot_Tradeoff

### ** Examples

## No test: 
# A typical example
set.seed(1234)
output = Plot_Tradeoff(
  data = TBILL,
  f1   = "SPFfor_Step1",
  f2   = "NCfor_Step1",
  y    = "Realiz1",
  m_set = seq(from = 1, to = 70, by = 10)
)
output[[1]] # The first element is a ggplot2 object of the size-power tradeoff.
output[[2]] # The second element is the underlying data used to construct the plot in element 1.

# An example with a user supplied loss function
# To use the mean absolute error as a loss function rather than a quadratic loss function
set.seed(1234)
output = Plot_Tradeoff(
  data = TBILL,
  f1   = "SPFfor_Step1",
  f2   = "NCfor_Step1",
  y    = "Realiz1",
  loss_function = function(f,y){ abs(f-y) },
  m_set = seq(from = 1, to = 50, by = 10)
)

# An example without (f1, f2, y). The function will take the first three columns and use them
set.seed(1234)
tmpdata = TBILL[, c("SPFfor_Step1", "NCfor_Step1", "Realiz1")] # data with [f1, f2, y]
Plot_Tradeoff(
  data = tmpdata,
  m_set = seq(from = 1, to = 50, by = 10)
)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Plot_Tradeoff", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
