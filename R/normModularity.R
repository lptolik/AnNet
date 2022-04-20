##---------------------------------------------------------------------
## Ref: Takemoto, K. & Kihara, K. Modular organization of cancer signaling networks is associated with patient survivability. Biosystems 113, 149–154 (2013).
##
## To allow the comparison of network modularity with networks of different
## size and connectivity.
##
## Used the normalized network modularity value
## Qm based on the previous studies by Parter et al., 2007, Takemoto, 2012,
## Takemoto, 2013, Takemoto and Borjigin, 2011, which was defined as:
## Qm = (Qreal-Qrand)/(Qmax-Qrand)
## Where Qreal is the network modularity of a real-world signaling network and,
## Qrand is the average network modularity value obtained from 10,000
## randomized networks constructed from its real-world network. Qmax was
## estimated as: 1 − 1/M, where M is the number of modules in the real network.
##
## Randomized networks were generated from a real-world network using the
## edge-rewiring algorithm (Maslov and Sneppen, 2002).
##---------------------------------------------------------------------



#' Calculates the normalized network modularity value.
#'
#' @param gg
#' @param alg
#' @param seed
#' @param Nint
#'
#' @return
#' @export
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' nm<-normModularity(gg,alg='louvain')
normModularity<-function(gg,alg=c('lec','wt','fc','infomap','louvain','sgG1','sgG2','sgG5'),seed=NULL,Nint=1000){
  cl<-getClustering(gg,alg)
  Qobs <- max(cl$modularity)

  ##--- get max number of communities
  Kmax <- max(as.numeric(cl$membership))

  ##--- define Qmax
  Qmax <- 1 - 1/Kmax

  ##--- generate random graph from our graph using rewring function,
  ##    preserving original graphs degree distribution.
  ##    Using 1000 rewing studies to get random modularity for graph with cl clustering
  Qrnd <- 0
  if(!is.null(seed)){
  set.seed(seed)
  }
  for( i in 1:Nint ){
    gg.rnd      = igraph::rewire(graph=gg,with=keeping_degseq(loops=FALSE,niter=100))
    cl.rnd = getClustering(gg.rnd,alg)
    Qrnd        = Qrnd + max(as.numeric(cl.rnd$modularity))
    rm(gg.rnd,cl.rnd)
  }

  ##--- random modularity for graph given cl clustering
  Qrnd = Qrnd/Nint

  ##--- normalised modularity for graph given cl clustering
  Qnorm = (Qobs-Qrnd)/(Qmax-Qrnd)
 return(Qnorm)
}
