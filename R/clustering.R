
cluster<-function(gg,alg=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5')){
  alg <- match.arg(alg)
  ids <- V(gg)$name
  lec<-function(gg){
    lec     <- igraph::leading.eigenvector.community(gg)
    ll      <- igraph::leading.eigenvector.community(gg, start=membership(lec))
  }
  cl<-switch(alg,
             lec=lec(gg),
             wt=igraph::walktrap.community(gg),
             fc=igraph::fastgreedy.community(gg),
             infomap=igraph::cluster_infomap(gg),
             louvain=igraph::cluster_louvain(gg),
             sgG1=igraph::spinglass.community(gg, spins=as.numeric(500),gamma=1),
             sgG2=igraph::spinglass.community(gg, spins=as.numeric(500),gamma=2),
             sgG5=igraph::spinglass.community(gg, spins=as.numeric(500),gamma=5)
               )
  cc       <- as.data.frame(names=cl$names,membership=cl$membership)
  attrib
  return(cc)
}

