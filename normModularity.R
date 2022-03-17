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


##--- declare clustering algorithms used 
##--- Clustering algorithms used 
alg <- "louvain"

##--- load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("PPI_Presynaptic.gml",format="gml")

##--- get node community assignment using louvain algorithm
louvain <- igraph::cluster_louvain(gg)

##--- get modularity value from running louvain algorithm
Qobs <- max(louvain$modularity)
                         
##--- get max number of communities
Kmax <- max(as.numeric(louvain$membership))

##--- define Qmax
Qmax <- 1 - 1/Kmax

##--- generate random graph from our graph using rewring function,
##    preserving original graphs degree distribution.
##    Using 1000 rewing studies to get random modularity for graph with louvain clustering
Qrnd <- 0
set.seed(1234)
Nints = 1000
for( i in 1:Nints ){
 gg.rnd      = igraph::rewire(graph=gg,with=keeping_degseq(loops=FALSE,niter=100))
 louvain.rnd = igraph::cluster_louvain(gg.rnd)
 Qrnd        = Qrnd + max(as.numeric(louvain.rnd$modularity))
 rm(gg.rnd,louvain.rnd)
}

##--- random modularity for graph given louvain clustering
Qrnd = Qrnd/Nints

##--- normalised modularity for graph given louvain clustering
Qnorm = (Qobs-Qrnd)/(Qmax-Qrnd)
