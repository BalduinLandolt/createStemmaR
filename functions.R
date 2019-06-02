


#
# Functions
# =========
#
# Here, the functions are declared, that are used in `script.R`.
#
#




#
# make_trees
# ----------
#
# This Function calculates different trees from a data set, using different algorithms
# 
# Parameters:
#   data (data.frame): an alignement of taxa
#   name (character): whatever the three is supposed to be called
#
# Returns:
#   list: a list of different tree-shaped interpretations of the input data
#

make_trees = function(data, name){
  
  # initialize the list eventually to be returned
  # add name to the list
  res = list(name = name)
  
  # generate paths for saving nexus and tree file to disk
  # (to be able to check results after run time) 
  nexus_path = paste("data/", name, ".nex", sep = "")
  tree_path = paste("data/", name, "_tree.tre", sep = "")
  
  # save nexus to disk
  save_nexus(data, f = nexus_path)
  
  # read nexus data from disk in propper nexus format
  d = read.nexus.data(nexus_path)
  res$nexus = d
  # the above is technically a workaround I came up with, to avoid data format problems.
  # storing the nexus file is just a nice side effect.
  
  
  # generic representation of data
  align = as.alignment(d)
  res$alignment = align
  res$character = as.character(res$alignment)
  
  
  # create phyDat object from nexus data
  phy = phyDat(d, type = "USER", levels = c("0","1","2","3","4","5","6","7","8","9","-"), ambiguity = c("?"))
  # store in result list
  res$phyDat = phy
  
  # Should be identical to phyDat
  # (so this is basically useless - but that way, the whole nexus-thing could be avoided.)
  phy2 = as.phyDat(align, type = "USER", levels = c("0","1","2","3","4","5","6","7","8","9","-"), ambiguity = c("?"))
  res$phy2 = phy2
  
  # calculate neighbour join tree from data
  tree = dist.ml(phy)
  nj_data = NJ(tree)
  # store in result list
  res$NJ = nj_data
  
  
  # NB for Max. Likelyhood, Max. Parsimony and Consensus Net:
  # Bootstrap Values of 20 (`bs = 20`) are too low for reliable results.
  # But having reasonable values (e.g. `100`) slows down the program a lot.
  
  # Maximum Likelyhood
  fit = pml(nj_data, phy)
  fit = optim.pml(fit, rearrangement="NNI")
  bs = bootstrap.pml(fit, bs=20, optNni=TRUE)
  res$maxLikely = bs
  
  # Maximum parsimony
  treeMP = pratchet(phy)
  treeMP = acctran(treeMP, phy)
  BStrees = bootstrap.phyDat(phy, pratchet, bs = 20)
  
  res$treeMP = treeMP
  res$maxParsimonyTrees = BStrees
  
  #consensus net
  cnt = consensusNet(bootstrap.phyDat(phy, FUN = function(x)nj(dist.hamming(x)), bs=20))
  res$consensusNet = cnt
  
  # neighbour net
  nnt = neighborNet(dist.hamming(phy))
  res$neighbourNet = nnt
  
  
  # save tree to file (reproducability)
  write.tree(nj_data, file=tree_path)
  
  return(res)
}


#
# plot_trees
# ----------
#
# This function plots the tree data to dendograms
#
# Parameters:
#   trees (list): a list of different tree data, as produced by `make_trees()`
#
# Returns:
#   NULL
#

plot_trees = function(trees){
  
  # plot different interpretations of the data
  
  plot(trees$NJ, main = "Neighbour Join - Type: Phylogram", sub = trees$name)
  plot(trees$NJ, main = "Neighbour Join - Type: Unrooted", type = "unrooted", sub = trees$name)
  plot(trees$NJ, main = "Neighbour Join - Type: Radial", type = "radial", sub = trees$name)
  
  plotBS(trees$NJ, trees$maxLikely, "unrooted", main="Maximum Likelyhood: Unrooted", sub = trees$name)
  plotBS(trees$NJ, trees$maxLikely, "phylogram", main="Maximum Likelyhood: Phylogram", sub = trees$name)
  
  plotBS(trees$treeMP, trees$maxParsimonyTrees, "unrooted", main="Maximum Parsimony: Unrooted", sub = trees$name)
  plotBS(trees$treeMP, trees$maxParsimonyTrees, "phylogram", main="Maximum Parsimony: Phylogram", sub = trees$name)
  
  plot(trees$consensusNet, "2D")
  title(main = "Consensus Net\n(2D Rendering)", sub = trees$name)
  
  plot(trees$neighbourNet, "2D")
  title(main = "Neighbour Net", sub = trees$name)
}



#
# save_nexus
# ----------
#
# This function saves given data to a given file in valid .nexus format
#
# Parameters:
#   data (data.frame): an alignment of taxa
#   f (character): path to the file, where data is supposed to be stored. (should end on `.nex`)
#
# Returns:
#   NULL
#
save_nexus = function(data, f){
  header = generate_header(data)
  body = generate_body(data)
  tail = get_tail()
  
  # save nexus file
  cat(header, body, tail, file=f, sep = "\n")
}




#
# generate_header
# ---------------
#
# This function generates the opening sequence of a nexus file.
# (Should only be called from `save_nexus()`)
# NB: The symboles are hard-coded here:
# 0:9 for regular symbols; - for gap; ? for missing.
#
# Parameters:
#   data (data.frame): an alignment of taxa; only used for the dimensions of the matrix.
#
# Returns:
#   character: said nexus header.
#
generate_header = function(data){
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
  header
}





#
# generate_body
# -------------
#
# This function generates the matrix sequence of a nexus file.
# (Should only be called from `save_nexus()`)
#
# Parameters:
#   data (data.frame): an alignment of taxa; from this data, the matrix of the nexus file is generated.
#
# Returns:
#   character: said nexus body
#
generate_body = function(data){
  data[[1]] = make_unique_names(data[[1]])
  lines = c()
  for (r in 1:nrow(data)){
    line = normalize_name(data[r,1])
    for (c in 2:ncol(data)){
      line = paste(line, as.character(data[r,c]), sep="")
    }
    lines = append(lines, line)
  }
  body = paste(lines, sep = "\n")
}




#
# make_unique_names
# -----------------
#
# This function makes sure that every name is only used once.
# If names occur multiple times, the last couple of characters get replaced by underscore plus an incrementing number.
#
# Parameters:
#   names (vector): names, potentially not unique.
#
# Returns:
#   vector: a vector of same length and order, but with unique names in it
#
make_unique_names = function(names){
  res = c()
  for (name in names) {
    name = normalize_name(name)
    i = 1
    while (name %in% res) {
      name = trimws(name)
      replacement = as.character(i)
      end = nchar(name)-nchar(replacement)-1
      name = substr(name, 1, end)
      name = paste(name, replacement, sep = "_")
      name = normalize_name(name)
      i = i+1
    }
    res = c(res, name)
  }
  res
}





#
# get_tail
# --------
#
# This function returns the end sequence of a nexus file.
# (Should only be called from `save_nexus()`)
#
# Parameters:
#   none
#
# Returns:
#   character: said nexus tail
#
get_tail = function(){
  tail = "; [...und hier endet es]
END; [beendet den Data-Block]"
  tail
}







#
# normalize_name
# --------------
#
# This function returns normalizes "species" names so that they are valid for a nexus file.
#
# A valid name is exactly 10 characters long, followed by a whitespace.
# The input gets shortened to 10 characters, if it's longer than that;
# if it's shorter, underscores are appended to make it 10 characters long.
#
# Blanks, hyphens and dots are replaced by Underscore.
# All other characters that don't match [A-Za-z0-9_] are removed without substitution.
# 
#
# Parameters:
#   x (character): input name to be normalized
#
# Returns:
#   character: normalization of the imput
# 
normalize_name = function(x){
  n = as.character(x)
  n = gsub(" ", "_", n)
  n = gsub("-", "_", n)
  n = gsub("\\.", "_", n)
  n = gsub("[^A-Za-z0-9_]", "", n)
  n = paste(n, "____________", sep = "")
  n = substr(n, start = 1, stop = 10)
  n = paste(n, " ", sep = "")
  n
}









#
# get_random_data
# ---------------
#
# This function generates a random data sample:
#    5-15 taxa (rows)
#    20-200 sequence length (columns)
#    in this alignment matrix, every data point is random.
#
# Of course, this data does not allow for sensible dendrograms; but it proofs that the script works.
# 
#
# Parameters:
#   None
#
# Returns:
#   data.frame: random data
# 
get_random_data = function(){
  # random dimensions
  rows = floor(runif(1, min=5, max=15))
  cols = floor(runif(1, min=20, max=200))
  
  frame = data.frame()
  
  # add names (first column)
  #nn = c()
  #for (r in 1:rows) {
  #  nn = append(nn, as.character(r))
  #}
  #frame$name = nn
  
  # add random data by column
  for (col in 1:cols){
    taxon = c()
    # do for every row in current column
    for (r in 1:rows) {
      i = floor(runif(1, min=0, max=9))
      frame[r,col] = as.character(i)
    }
  }
  
  frame = cbind(rownames(frame), frame)
  
  # return data frame
  frame
}




make_tree = function(d, name){
  
  nexus_path = paste("data/", name, ".nex", sep = "")
  tree_path = paste("data/", name, "_tree.tre", sep = "")
  
  save_nexus(data, f = nexus_path)
  
  
  # read nexus data
  d = read.nexus.data(nexus_path)
  
  # create phyDat object from nexus
  phy = phyDat(d, type = "USER", levels = c("?","-","0","1","2","3","4","5","6","7","8","9"))
  
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
  
  
  
  # save tree to file
  write.tree(nj_data, file=tree_path)
  
  
  plot(nj_data, main = "Plot with Title 1")
  plot(nj_data, main = "Plot with Title 2", type = "unrooted")
}

