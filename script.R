
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

source("utils_phylo.R")
source("tree_maker.R")



#
# Make Tree with Function
# =======================
#
#


data = read.csv("data/trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)
make_tree(d, name = "trivial")
data = get_random_data()
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




