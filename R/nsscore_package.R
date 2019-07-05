#' nsscore: Normal score transformation and back transformation for kriging
#'
#' \pkg{nsscore} Experimental R functions for performing normal score transforms and back transforms
# on data. This is useful when conducting geostatistical Gaussian simulation.
#
# Kriging with normal scores may be necessary though not sufficient for ensuring that a problem
# reflects the multiGaussian assumptions of sequential Gaussian simulation. Checks are available;
# see Goovaerts (1999) chapter 7.
#
# Functions were developed that are loosely based on fortran routines in GSLIB v.2; see Deutsch &
# Journel, 1998. Errors, bugs, and deviations are all due to my coding and implementation decisions,
# however.
#
# The back transform implemented here interpolates linearly between data. Extrapolation is also
# linear. See backtr() code for details and options. This back transform implementation is
# problematic for two reasons:
# 1. How to evaluate ties (it's a rank transform)
# 2. Extrapolating for small and large values.
#
# Valuable discussions on AI-Geostats archives as well as in the formal literature consider these
# issues, and possibly much better ideas are out there (e.g. Saito & Goovaerts (2000))
#
# That said, extrapolation decisions in this code are largely theory-free and the user is warned to
# treat such values with suspicion.
#
# These functions are experimental and are not written to be robust. In particular, function inputs
# are not tested for validity, there's no error handling, etc. etc. It is provided in hopes that it
# will be stimulating.
#
# For examples, try:
# source('nscore.R')   # loads the functions, runs nothing
# example1.nscore()
# example2.nscore()
#
# References
# Deutsch, C.V. and Journel, A.G. (1998) GSLIB: Geostatistical
# Software Library and User's Guide. New York:Oxford.
#
# Goovaerts, P. (1997) Geostatistics for Natural Resources
# Evaluation. New York:Oxford.
#
# Dubois, G. (2001) AI-GEOSTATS: SUMMARY: Nscore transform &
# kriging of log normal data sets.
# http://www.mail-archive.com/ai-geostats@jrc.it/msg00152.html
#
# Saito, H. & Goovaerts, P. (2000) Geostatistical interpolation of
# positively skewed and censored data in a dioxin-contaminated site.
# Environmental Science & Technology 34(19): 4228-4235.
#
#' @author Ashton Shortridge, May/June, 2008
#' @docType package
#' @name nsscore
#' @examples
#' # Example 1 ----
#' # This function illustrates the use of nscore and backtr on a nonspatial synthetic example.
#' rain <- runif(1000, 0, 7)
#' rain <- (rain)^2
#' rain[rain < 0] <- 0 # A rainfall-like distribution - right skewed
#' hist(rain, main = "Hypothetical population of precipitation estimates") # Yep, right skewed.
#' rsam <- sample(rain, 100) # Extract a sample from rain
#' hist(rsam) # Eerily similar, yet different
#'
#' rsam.sc <- nscore(rsam) # normalize the sample
#' plot(rsam.sc$trn.table$nscore, rsam.sc$trn.table$x, main = "Sample rainfall normal scores")
#' lines(approx(rsam.sc$trn.table$nscore, rsam.sc$trn.table$x)) # approx is a linear interp
#'
#' # Back transform the nscored transform
#' rsam.bak <- backtr(rsam.sc$nscore, rsam.sc, tails = "separate")
#' cor(rsam, rsam.bak) # Nice.
#'
#' # Nscore and backtransform the rain population data using the sample transform table
#' rain.sc <- nscore(rain)
#' rain.back <- backtr(rain.sc$nscore, rsam.sc, tails = "none")
#'
#' cor(rain, rain.back) # about 0.99; pretty close. Note the issues in the tails.
#' summary(cbind(rain, rsam, rain.back))
#' par(mfrow = c(2, 2))
#' hist(rain, main = "Hypothetical population of precipitation estimates")
#' hist(rain.sc$nscore, main = "Population scores")
#' hist(rain.back, main = "Back-transformed precip estimates using sample transform table")
#' qqplot(rain, rain.back, cex = 0.6, main = "Distribution Comparison")
#'
#' # Now, can we draw 1000 standard normal values and use rsam.sc to reproduce rain?!
#' rain.sim.sc <- rnorm(1000, 0, 1)
#' rain.sim <- backtr(rain.sim.sc, rsam.sc, tails = "equal") # use the sample backtransform table
#' hist(rain, main = "Hypothetical population of precipitation estimates")
#' hist(rain, main = "Back-transformed simulated normal scores using the sample transform table ")
#' qqplot(rain, rain.sim, cex = 0.6, main = "Distribution Comparison")
#' print(summary(cbind(rain, rsam, rain.sim))) # summary stats....
#' par(mfrow = c(1, 1))
#'
#' # Example 2 ----
#' # An example using gstat library simulation routines and the Meuse dataset.
#' library(gstat)
#' library(sp)
#' data(meuse)
#' coordinates(meuse) <- ~ x + y
#' hist(log(meuse$zinc), main = "ln(Zinc): Not too Gaussian")
#' v <- variogram(log(zinc) ~ 1, meuse)
#' m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
#' print(plot(v, model = m, main = "Variogram model on ln(Zn)"))
#' set.seed(131)
#' data(meuse.grid)
#' gridded(meuse.grid) <- ~ x + y
#' sim.log <- krige(
#'   formula = log(zinc) ~ 1, meuse, meuse.grid, model = m,
#'   nmax = 15, beta = 5.9, nsim = 4
#' )
#' # show all 4 simulations
#' print(spplot(sim.log, main = "Conditional SK Simulation on ln(Zn)"))
#'
#' # Back transform
#' sim <- sim.log
#' sim$sim1 <- exp(sim$sim1)
#' sim$sim2 <- exp(sim$sim2)
#' sim$sim3 <- exp(sim$sim3)
#' sim$sim4 <- exp(sim$sim4)
#' print(spplot(sim, main = "Back-transformed ln(Zn) simulations"))
#'
#' # Try a normal score transform instead
#' meuse.zinc.sc <- nscore(meuse$zinc)
#' meuse$zinc.sc <- meuse.zinc.sc$nscore
#' hist(meuse$zinc.sc, main = "Zinc normal score")
#' v <- variogram(zinc.sc ~ 1, meuse)
#' m <- fit.variogram(v, vgm(1, "Sph", 300, 1))
#' print(plot(v, model = m, main = "Variogram model on normal scores Zn"))
#' set.seed(131)
#' sim.ns <- krige(
#'   formula = zinc.sc ~ 1, meuse, meuse.grid, model = m,
#'   nmax = 15, beta = 0, nsim = 4
#' )
#' # show all 4 ns simulations
#' print(spplot(sim.ns, main = "Conditional SK Simulation on Zn scores"))
#' hist(sim.ns$sim1, main = "Histogram of simulation 1 - pretty normal")
#' sim.ns$sim1 <- backtr(sim.ns$sim1, meuse.zinc.sc, tails = "separate")
#' sim.ns$sim2 <- backtr(sim.ns$sim2, meuse.zinc.sc, tails = "separate")
#' sim.ns$sim3 <- backtr(sim.ns$sim3, meuse.zinc.sc, tails = "separate")
#' sim.ns$sim4 <- backtr(sim.ns$sim4, meuse.zinc.sc, tails = "separate")
#' print(spplot(sim.ns, main = "Back-transformed Zn nscore simulations"))
#' print("Original Zn data")
#' print(summary(meuse$zinc))
#' print("Back transformed Zn normal score simulation")
#' print(summary(sim.ns$sim1))
#' print("Back transformed ln(Zn) simulation")
#' print(summary(sim$sim1))
NULL
