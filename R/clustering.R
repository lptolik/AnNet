
#' Calculate cluster memberships for the graph.
#'
#' @param gg
#' @param alg
#'
#' @return
#' @export
#'
#' @examples
calcMembership<-function(gg,alg=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5','spectral')){
  ids <- V(gg)$name
  cl<-getClustering(gg,alg)
  if(!is.null(cl)){
    cc       <- data.frame(names=cl$names,membership=cl$membership)
  }else{
    cc <- data.frame(names='names',membership=0)[FALSE,]
  }
  return(cc)
}

#' Calculate memberships for all clustering algorithms and store them on the
#' graph vertices.
#'
#' @param gg
#'
#' @return
#' @export
#'
#' @examples
calcAllClustering<-function(gg){
  ids <- V(gg)$name
  cnames<-c('ID','lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5','spectral')
  m      <- matrix(NA, ncol=length(cnames), nrow=length(ids))
  colnames(m)<-cnames
  m[,1]<-ids
  for(ai in 2:length(cnames)){
    an<-colnames(m)[ai]
    cm<-calcMembership(gg,an)
    if(dim(cm)[1]>0){
      m[,ai]<-as.character(cm$membership)
      mod<-modularity(gg,cm$membership)
      gg<-set.graph.attribute(gg,an,mod)
    }
  }
  ggm<-applpMatrixToGraph(gg,m)
  return(ggm)
}

#' Calculate memberships for particular clustering algorithms and store them on the
#' graph vertices.
#'
#'
#' @param gg
#' @param alg algorithm to apply
#'
#' @return
#' @export
#'
#' @examples
calcClustering<-function(gg,alg){
  cl<-getClustering(gg,alg)
  if(!is.null(cl)){
    ids <- V(gg)$name
    m      <- matrix(NA, ncol=2, nrow=length(ids))
    colnames(m)<-c('ID',alg)
    m[,1]<-ids
    m[,2]<-as.character(cl$membership)
    ggm<-applpMatrixToGraph(gg,m)
    mod<-modularity(ggm,cl$membership)
    ggm<-set.graph.attribute(ggm,alg,mod)
    return(ggm)
  }else{
    return(gg)
  }
}

#' Get clustering results for the graph.
#'
#' @param gg
#' @param alg
#'
#' @return
#' @export
#'
#' @examples
getClustering<-function(gg,alg=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5','spectral')){
  alg <- match.arg(alg)
  lec<-function(gg){
    lec     <- igraph::leading.eigenvector.community(gg)
    ll      <- igraph::leading.eigenvector.community(gg, start=membership(lec))
  }
  cl<-try(switch(alg,
                 lec=lec(gg),
                 wt=igraph::walktrap.community(gg),
                 fc=igraph::fastgreedy.community(gg),
                 infomap=igraph::cluster_infomap(gg),
                 louvain=igraph::cluster_louvain(gg),
                 sgG1=igraph::spinglass.community(gg,
                                                  spins=as.numeric(500),gamma=1),
                 sgG2=igraph::spinglass.community(gg,
                                                  spins=as.numeric(500),gamma=2),
                 sgG5=igraph::spinglass.community(gg,
                                                  spins=as.numeric(500),gamma=5),
                 spectral=rSpectral::spectral_igraph_communities(gg)
  ))
  if(inherits(cl, "try-error")){
    warning('Clustering calculations for algorithm "',alg,
            '" failed. NULL is returned')
    return(NULL)
  }
  return(cl)
}

#' Matrix of cluster characteristics
#'
#' @param gg graph to analyse
#' @param att vector of attribute names that contains membership data
#'
#' @return matrix of clustering characteristics
#' @export
#'
#' @examples
#'
clusteringSummary<-function(gg,att=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5','spectral')){
  attN<-vertex_attr_names(gg)
  idx<-match(attN,att)
  clusterings<-attN[!is.na(idx)]
  res<-list()
  for(c in clusterings){
    cmem<-as.numeric(vertex_attr(gg,c))
    mod<-modularity(gg,cmem)
    Cn<-table(cmem)
    C <- length(Cn)
    Cn1<-length(which(Cn==1))
    Cn100<-length(which(Cn>=100))
    summary(as.vector(Cn))->s
    names(s)<-paste(names(s),'C')
    sgraphs<-lapply(names(Cn),getClusterSubgraphByID,gg=gg,mem=cmem)
    ug <- disjoint_union(sgraphs)
    mu<- 1-ecount(ug)/ecount(gg)
    r1<-c(mod,C,Cn1,Cn100,mu)
    names(r1)<-c('mod','C','Cn1','Cn100','mu')
    res[[c]]<-c(r1,s)
  }
  return(do.call(rbind,res))
}
