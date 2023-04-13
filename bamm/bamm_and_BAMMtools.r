#' To examine diversification rate shift based on birth-death model using BAMM and BAMMtools
#' Generated 2023/03/21 by Shusuke Kawakubo
#' Last updated 2023/04/13 by Shusuke Kawakubo
#'
#' FYI please refer to the following site: http://bamm-project.org/index.html

# load packages
library("BAMMtools")
library("phytools")
library("coda")
library("magrittr")
library("RColorBrewer")

# First make sure to export beast output mcc tree as newick format through FigTree
# read tree
P1_mcc_total <- read.tree("./MCC_tree/P1_TotalIsolates_MCC_newick.tre")

# make sure to fix the tree ultrametric
P1_mcc_total <- ladderize(force.ultrametric(P1_mcc_total, method = "extend"), right = FALSE)
is.ultrametric(P1_mcc_total)
plot(P1_mcc_total)
write.tree(P1_mcc_total, file = "./bamm_output/P1_total/P1_mcc_total_ultrametric.tre")

# set the appropriate prior distribution, assuming that current data set contain 1% of global population
setBAMMpriors(P1_mcc_total, total.taxa = 17900, traits = NULL,
              outfile = "./bamm_output/P1_total/P1_total_myPriors.txt")

# generate the control file for bamm analysis with inferred prior based on the "setBAMMpriors" function
generateControlFile(file = "./bamm_output/P1_total/P1_total_controlFile.txt", type = "diversification", params = list(
  treefile = 'P1_mcc_total_ultrametric.tre',
  globalSamplingFraction = '0.01',
  numberOfGenerations = '5000000',
  mcmcWriteFreq = "5000",
  eventDataWriteFreq = "5000",
  printFreq = "5000",
  overwrite = '1',
  lambdaInitPrior = '19.1906569899343',
  lambdaShiftPrior = '0.00140873513838022',
  muInitPrior = '19.1906569899343',
  expectedNumberOfShifts = '1'))

# open terminal and run a following command in the directry
# $ bamm -c ./bamm_output/P1_total/P1_total_myPriors.txt
# Optionally, reloading ".zprofile" may be required if you are using terminal imbedded in Rstudio. If so, run the following:$ source .zprofile


# import BAMM output as "bammdata" object
P1_total_edata <- getEventData(P1_mcc_total, eventdata = "./bamm_output/P1_total/event_data.txt", burnin=0.1)

# assess the MCMC convergence
mcmcout <- read.csv("./bamm_output/P1_total/mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

# discard burn-in states as appropriate proportion
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# check effective sample size (ESS) values if >200 is achieved
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)


# plot phylorate
summary(P1_total_edata)
plot.bammdata(P1_total_edata, lwd=1, legend=T)
title("P1_total")

# to check the specific MCMC sample, for example, 25th sample from the whale posterior:
# index <- 25
# e2 <- subsetEventData(P1_total_edata, index = index)
# plot.bammdata(e2, lwd=2)
# addBAMMshifts(e2, cex=2)

# Bayesian credible sets of shift configurations
css <- credibleShiftSet(P1_total_edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
plot.credibleshiftset(css, lwd = 1, legend = T)

best <- getBestShiftConfiguration(P1_total_edata, expectedNumberOfShifts=1)
plot.bammdata(best, lwd = 1, legend=T)
addBAMMshifts(best, cex=2.5)
title("P1_total")

# plot phylorate
summary(P1_total_edata)
plot.bammdata(P1_total_edata, lwd=1, legend=T)
addBAMMshifts(best, cex=2.5)
title("P1_total")

# now clade specific rate
# you can identify node numbers using plot.phylo and nodelabels from the ape package
# First find the node number of your interest
plot.phylo(P1_mcc_total, cex = .5)
nodelabels(cex = .5)

# Clade-specific evolutionary rates
# First plot all rates
allrates <- getCladeRates(P1_total_edata)

# plot birth rate
#hist(allrates$lambda, breaks=12, col="purple")
mean(allrates$lambda)
#median(allrates$lambda)
#quantile(allrates$lambda, c(0.05, 0.95))
d_P1_total_lambda <- density(allrates$lambda)
plot(d_P1_total_lambda, main="P1 total birth rate", xlim = c(0.05,0.4), ylim = c(0, 20))
polygon(d_P1_total_lambda, col="#3a257f7F", border= "#3a257f") # "4C" = 30%, "7F" = 50%, "99" = 60%, "B2" = 70% transparency

# check leek clade
P1_total_leek <- getCladeRates(P1_total_edata, node=214)

d_P1_total_leek_lambda <- density(P1_total_leek$lambda)
#plot(d_P1_total_leek, main="P1 total leek birth rate", xlim = c(0.05,0.4), ylim = c(0, 20))
polygon(d_P1_total_leek_lambda, col="#a7c73f7F", border="#a7c73f")
text(x = .33, y = 15, cex = 1.5, labels = "mean birth rate")
text(x = .33, y = 13, cex = 1.5, labels = expression(paste(lambda [all] == "0.200")))
text(x = .33, y = 11, cex = 1.5, labels = expression(paste(lambda [leek] == "0.205")))
#hist(P1_total_leek$lambda)
mean(P1_total_leek$lambda)
#median(P1_total_leek$lambda)
#quantile(P1_total_leek$lambda)

# next plot death rate
#hist(allrates$lambda, breaks=12, col="purple")
mean(allrates$mu)
#median(allrates$lambda)
#quantile(allrates$lambda, c(0.05, 0.95))
d_P1_total_mu <- density(allrates$mu)
plot(d_P1_total_mu, main="P1 total death rate", xlim = c(0.05,0.4), ylim = c(0, 20))
polygon(d_P1_total_mu, col="#3a257f7F", border= "#3a257f") # "4C" = 30%, "7F" = 50%, "99" = 60%, "B2" = 70% transparency

# check leek clade
P1_total_leek <- getCladeRates(P1_total_edata, node=214)

d_P1_total_leek_mu <- density(P1_total_leek$mu)
#plot(d_P1_total_leek, main="P1 total leek birth rate", xlim = c(0.05,0.4), ylim = c(0, 20))
polygon(d_P1_total_leek_mu, col="#a7c73f7F", border="#a7c73f")
text(x = .33, y = 15, cex = 1.5, labels = "mean death rate")
text(x = .33, y = 13, cex = 1.5, labels = expression(paste(mu [all] == "0.179")))
text(x = .33, y = 11, cex = 1.5, labels = expression(paste(mu [leek] == "0.183")))
#hist(P1_total_leek$lambda)
mean(P1_total_leek$mu)
#median(P1_total_leek$lambda)
#quantile(P1_total_leek$lambda)

# checl S-type clade
P1_total_Stype <- getCladeRates(P1_total_edata, node=286)

d_P1_total_Stype <- density(P1_total_Stype$lambda)
#plot(d_P1_total_leek, main="P1 total leek birth rate", xlim = c(0.05,0.4), ylim = c(0, 20))
polygon(d_P1_total_Stype, col=rgb(0, 0, 0, alpha=0.1), border="pink")

#' optional analysis
#'
#' plotRateThroughTime(P1_total_edata, ratetype="speciation")
#' cmat <- getCohortMatrix(P1_total_edata)
#' cohorts(cmat, P1_total_edata)
#' cst <- cumulativeShiftProbsTree(P1_total_edata)
#' plot.phylo(cst)
#'

