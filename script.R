
# Run instructions:
#
# - set working directory to where the data lies, or adjust paths accordingly.


# imports
#
# Needs packages:
# - anchors
# -

# for easy replacement of values in data frames
library(anchors)


# load data
data = read.csv("trivial.csv", header = T, sep = ";", stringsAsFactors = F)
head(data)


data=replace.value(data, names = colnames(data) ,from = "-", to = "-1")
data=replace.value(data, names = colnames(data) ,from = "?", to = "-99")
head(data)
str(data)

#data[] = tapply(data, function(x) as.numeric(as.character(x)))
names = data[,1]
data[1] = NULL
data = as.data.frame(sapply(data, as.numeric))
row.names(data) = names
#data = cbind(names, data)

head(data)
str(data)
data





# bullshit so far...


data = read.csv("trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)
head(data)

# collate data to sequence
collate = function(x){
  d_v = unlist(x)
  d_v = sapply(d_v, as.character)
  # normalize text name
  res = paste(d_v, sep = "")
  res
}

header = paste(nrow(data), ncol(data), sep = " ")
#d = sapply(data, collate)
lines = c()
for (r in 1:nrow(data)){
  # TODO normalize text name
  #line = c(data[r,1])
  line = data[r,1]
  for (c in 2:ncol(data)){
    #line = append(line, as.character(data[r,c]))
    line = paste(line, as.character(data[r,c]), sep="")
  }
  lines = append(lines, line)
  
  
  #li = c(data[i,])
  #print(li)
  #lines = c(lines, li)
}
lines[1]
dd = c()
for (variable in vector) {
  
}
str(d)
summary(d)
