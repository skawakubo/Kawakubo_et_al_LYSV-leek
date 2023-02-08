#' To perform Discriminant analysis of principal components (DAPC) analysis
#' Generated 2022/11/08 by Shusuke Kawakubo
#' Last updated 2023/02/08 by Shusuke Kawakubo

#' reference
#' https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

# load packages
library("adegenet")
library("poppr")
library("seqinr")

# import data
Data<- read.alignment(file="seq.fasta", format="fasta", forceToLower = F)
population <- read.table(file="pop.txt",head=T)
Dataind <- alignment2genind(Data, pop=population$pop, exp.char=c("a","t","g","c"), na.char="-",polyThres=1/100)


# set color palette of your choice
grpName  <- c("garlic", "leek", "rakkyo","HybridGarlic","elephant")
myCol <- c("#a94293", "#a7c73f", "#9485b6", "#4a3a3a", "#efba99")


# perform an initial DAPC analysis##
dapcdata <- dapc(Dataind, var.contrib = TRUE, scale = FALSE, n.pca = NULL, n.da =NULL)
scatter(dapcdata, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

# cross validation to look for the range of principal component (PC)
set.seed(20221108)
pcx <- xvalDapc(tab(Dataind, NA.method = "mean"), pop(Dataind))
pcx[-1]

# find the best PC value; PC_best
set.seed(20221108)
	pcx <- xvalDapc(tab(Dataind, NA.method = "mean"), pop(Dataind),n.pca = 35:45, n.rep = 1000, parallel = "multicore", ncpus = 4L)
# A=PC_best-5，B=PC_best+5，i.e, if PC_best=10, then A=10-5=5, B=10+5=15, so change 'A:B' above to '5:15'
pcx[-1]

# plot
scatter(pcx$DAPC, col = myCol, cex = 2, legend = TRUE,clabel = FALSE, posi.leg = "topleft", scree.pca = TRUE, posi.pca = "topright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

scatter(pcx$DAPC,1,1, col = myCol, legend = TRUE,txt.leg=paste(grpName))

