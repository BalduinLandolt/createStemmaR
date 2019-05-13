
# Run instructions:
#
# - set working directory to where the data lies, or adjust paths accordingly.


# imports
#
# Needs packages:
# - anchors
# -

library(ape)





data = read.csv("trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)
head(data)


normalize_name = function(x){
  n = as.character(x)
  n = paste(n, "          ", sep = "")
  n = substr(n, start = 1, stop = 9)
  n = paste(n, " ", sep = "")
  n
}


header = paste(nrow(data), ncol(data)-1, sep = " ")
#d = sapply(data, collate)
lines = c()
for (r in 1:nrow(data)){
  line = normalize_name(data[r,1])
  for (c in 2:ncol(data)){
    line = paste(line, as.character(data[r,c]), sep="")
  }
  lines = append(lines, line)
}

# TODO make names unique


cat(header, lines, file="tax.txt", sep = "\n")

d = read.dna("tax.txt", format = "sequential")
str(d)
d
