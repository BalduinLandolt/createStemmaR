
# General Information:
# ====================
#
# For a more comprehensive introduction to this script, see 'Readme.md' on:
# https://github.com/BalduinLandolt/createStemmaR
#


# Run instructions:
# =================
#
# - set working directory to where the data lies, or adjust paths accordingly.
# - install required packages
# - run script step by step


# imports:
# ========
#
# Required packages:
# - ape
#   [run: install.packages("ape") ]
# - phangorn
#   [run: install.packages("phangorn") ]





# Script starts here
# ==================




# load packages
library(ape)
library(phangorn)

# load functions
source("functions.R")




# load pre-made trivial data sample
data = read.csv("data/trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)

# print head of data
head(data)

# from that data, calculate trees using several different algorythms,
# and plot the according dendrograms
# (method that does that is specified below)
make_tree(d, name = "trivial")





# generate a random sample of data
# (you'll get complete nonsense, of course.)
data = get_random_data()

# print head of data
head(data)

# calculate and plot the same dendrograms from this data
make_tree(data, name = "random")




#
# Make Tree of trivial Data Sample
# ================================
#
#



# load trivial sample data
data = read.csv("data/trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)
head(data)


safe_nexus(data, f = "data/tax.nex")


# read nexus data
d = read.nexus.data("data/tax.nex")
str(d)
d

# create phyDat object from nexus
phy = phyDat(d, type = "USER", levels = c("?","-","0","1","2","3","4","5"))
str(phy)
summary(phy)

# calculate trees from data
# -------------------------
# Distance Matrix
tree = dist.ml(phy)

# Neighbour Join
nj_data = NJ(tree)

# plot tree (in several different looks)
plot.phylo(nj_data, use.edge.length=FALSE, cex=0.75)
plot.phylo(nj_data, use.edge.length=TRUE)
plot.phylo(nj_data, type = "unrooted", lab4ut = "axial")
plot.phylo(nj_data, type = "unrooted")
plot.phylo(nj_data, type = "radial")


# max pars
pars_data = pratchet(phy)

# plot tree (in several different looks)
plot.phylo(pars_data, use.edge.length=FALSE, cex=0.75)
plot.phylo(pars_data, use.edge.length=TRUE)
plot.phylo(pars_data, type = "unrooted", lab4ut = "axial")
plot.phylo(pars_data, type = "unrooted")
plot.phylo(pars_data, type = "radial")


bt = bootstrap.phyDat(phy,FUN = function(x)nj(dist.hamming(x)), bs=50)
#consensus
cnt = consensusNet(bt)
plot(cnt)



# safe tree to files
write.tree(nj_data, file="data/trivial_tree.tre")

# -> old notes... might still be of use
# library(anchors)
# data=replace.value(data, names = colnames(data) ,from = "-", to = "-1")








#
# Make Tree of Random Data Sample
# ================================
#
#



# load trivial sample data
data = get_random_data()
head(data)


safe_nexus(data, f = "data/random.nex")


# read nexus data
d = read.nexus.data("data/random.nex")
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





# max pars
pars_data = pratchet(phy)

# plot tree (in several different looks)
plot.phylo(pars_data, use.edge.length=FALSE, cex=0.75)
plot.phylo(pars_data, use.edge.length=TRUE)
plot.phylo(pars_data, type = "unrooted", lab4ut = "axial")
plot.phylo(pars_data, type = "unrooted")
plot.phylo(pars_data, type = "radial")


bt = bootstrap.phyDat(phy,FUN = function(x)nj(dist.hamming(x)), bs=50)
#consensus
cnt = consensusNet(bt)
plot(cnt, "2D")



# safe tree to file
write.tree(nj_data, file="data/random_tree.tre")



