
# Run instructions:
#
# - set working directory to where the data lies, or adjust paths accordingly.


# imports
#
# Needs packages:
# - ape
#   [run: install.packages("ape") ]
# - phangorn
#   [run: install.packages("phangorn") ]


library(ape)
library(phangorn)




# load trivial sample data
data = read.csv("data/trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)
head(data)


source("utils_phylo.R")
safe_nexus(data, f = "data/tax.nex")


# read nexus data
d = read.nexus.data("data/tax.nex")
str(d)
d

# create phyDat object from nexus
phy = phyDat(d, type = "USER", levels = c("?","-","0","1","2","3","4","5"))
str(phy)
summary(phy)

# calculate tree from data
tree = dist.ml(phy)
nj_data = NJ(tree)

# plot tree (in several different looks)
plot.phylo(nj_data, use.edge.length=FALSE, cex=0.75)
plot.phylo(nj_data, use.edge.length=TRUE)
plot.phylo(nj_data, type = "unrooted", lab4ut = "axial")
plot.phylo(nj_data, type = "unrooted")
plot.phylo(nj_data, type = "radial")

# safe tree to file
write.tree(nj_data, file="data/tree.tre")




# -> old notes... might still be of use
# library(anchors)
# data=replace.value(data, names = colnames(data) ,from = "-", to = "-1")


