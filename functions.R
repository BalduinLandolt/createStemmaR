


#
# Functions
# =========
#
#
#
#



make_trees = function(data, name){
  
  res = list(name = name)
  
  
  nexus_path = paste("data/", name, ".nex", sep = "")
  tree_path = paste("data/", name, "_tree.tre", sep = "")
  
  safe_nexus(data, f = nexus_path)
  
  
  # read nexus data
  d = read.nexus.data(nexus_path)
  
  # create phyDat object from nexus
  phy = phyDat(d, type = "USER", levels = c("?","-","0","1","2","3","4","5","6","7","8","9"))
  res$phyDat = phy
  
  # calculate tree from data
  tree = dist.ml(phy)
  nj_data = NJ(tree)
  res$NJ = nj_data
  
  
  
  # max pars
  pars_data = pratchet(phy)
  res$maximum_parsimony = pars_data
  
  
  
  bt = bootstrap.phyDat(phy,FUN = function(x)nj(dist.hamming(x)), bs=50)
  res$bootstrap50 = bt
  #consensus
  cnt = consensusNet(bt)
  res$consensusNet = cnt
  
  
  
  # safe tree to file
  write.tree(nj_data, file=tree_path)
  
  return(res)
}



plot_trees = function(trees){
  
  #plot(trees$phyDat, main = "PhyDat")
  
  plot(trees$NJ, main = "Neighbour Join - Type: Phylogram", sub = trees$name)
  plot(trees$NJ, main = "Neighbour Join - Type: Unrooted", type = "unrooted", sub = trees$name)
  plot(trees$NJ, main = "Neighbour Join - Type: Radial", type = "radial", sub = trees$name)
  
  plot(trees$maximum_parsimony, main = "Maximum Parsimony - Type: Phylogram", sub = trees$name)
  plot(trees$maximum_parsimony, main = "Maximum Parsimony - Type: Unrooted", type = "unrooted", sub = trees$name)
  plot(trees$maximum_parsimony, main = "Maximum Parsimony - Type: Radial", type = "radial", sub = trees$name)
  
  #plot(trees$consensusNet, main = "Consensus Net", sub = trees$name)
  plot(trees$consensusNet, "2D", main = "Consensus Net", sub = trees$name)
  
}




# function to create nexus as string
safe_nexus = function(data, f){
  header = generate_header(data)
  body = generate_body(data)
  tail = get_tail()
  
  # safe nexus file
  cat(header, body, tail, file=f, sep = "\n")
}




# generate nexus file header
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




normalize_name = function(x){
  n = as.character(x)
  n = paste(n, "____________", sep = "")
  n = substr(n, start = 1, stop = 10)
  n = paste(n, " ", sep = "")
  n
}




# generate data alignment string from data
generate_body = function(data){
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




get_tail = function(){
  tail = "; [...und hier endet es]
END; [beendet den Data-Block]"
  tail
}



# TODO make names unique




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
  
  safe_nexus(data, f = nexus_path)
  
  
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
  
  
  
  # safe tree to file
  write.tree(nj_data, file=tree_path)
  
  
  plot(nj_data, main = "Plot with Title 1")
  plot(nj_data, main = "Plot with Title 2", type = "unrooted")
}

