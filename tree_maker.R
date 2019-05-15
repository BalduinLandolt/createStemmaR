
make_tree = function(d, name){
  
  nexus_path = paste("data/", name, ".nex", sep = "")
  tree_path = paste("data/", name, "_tree.tre", sep = "")
  
  safe_nexus(data, f = nexus_path)
  
  
  # read nexus data
  d = read.nexus.data(nexus_path)
  
  # create phyDat object from nexus
  phy = phyDat(d, type = "USER", levels = c("?","-","0","1","2","3","4","5"))
  
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
}