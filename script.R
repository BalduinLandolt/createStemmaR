
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


normalize_name = function(x){
  n = as.character(x)
  n = paste(n, "          ", sep = "")
  n = substr(n, start = 1, stop = 9)
  n = paste(n, " ", sep = "")
  n
}

# generate nexus file header
header1 = "#NEXUS
BEGIN data;[eröffnet den Data-Block]
Dimensions ntax="
header2 = " nchar="
header3="; [Definiert die Größe des Alignments]
Format
datatype=standard
missing=?
gap=- [Definiert den Datentyp und Symbole für fehlende Daten (?) und gaps (-)]
Symbols=\"0 1 2 3 4 5 6 7 8 9\";
Matrix [hier beginnt das Alignment...]"
header = paste(header1, nrow(data), header2, ncol(data)-1, header3, sep = "")


# generate data alignment string from data
lines = c()
for (r in 1:nrow(data)){
  line = normalize_name(data[r,1])
  for (c in 2:ncol(data)){
    line = paste(line, as.character(data[r,c]), sep="")
  }
  lines = append(lines, line)
}

tail = "; [...und hier endet es]
END; [beendet den Data-Block]"

# TODO make names unique

# safe nexus file
cat(header, lines, tail, file="data/tax.nex", sep = "\n")

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

