##---Semi-local Centrality (Cl)
##   Identifying influential nodes in complex networks, D. Chen et al., Physica A, 2012
Semilocal <- function(gg){

  N    <- length(V(gg)$name)
  meas <- matrix(0, nrow=N, ncol=3)

  for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- igraph::neighbors(gg,v=ids,mode="all")

    if( length(neig) > 0 ){

      for( w in 1:length(neig) ){
        neig <- c(neig,igraph::neighbors(gg,v=as.character(V(gg)$name[neig[w]]),mode="all"))
      }

      neig <- unique(neig)

      meas[i,1] <- length(neig)-1

    }

  }

  for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- igraph::neighbors(gg,v=ids,mode="all")

    meas[i,2] <- sum(as.numeric(meas[neig,1]))

  }


  for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- igraph::neighbors(gg,v=ids,mode="all")

    meas[i,3] <- sum(as.numeric(meas[neig,2]))

  }

  return(as.numeric(meas[,3]))

}

fSemilocal<-function(gg){
  N    <- length(V(gg)$name)
  meas <- matrix(0, nrow=N, ncol=3)
  meas[,1]<-ego_size(gg,order = 2,mode='all')-1
  neigSum<-function(i,graph,vec){
    neig <- igraph::neighbors(graph,v=i,mode="all")
    return(sum(vec[neig]))
  }
  meas[,2]<-sapply(1:N,neigSum,graph=gg,vec=meas[,1])
  meas[,3]<-sapply(1:N,neigSum,graph=gg,vec=meas[,2])
  return(as.numeric(meas[,3]))
}

##calculate the mean and sd of the shortest paths for each gene
calShorestPaths <- function(gg){

  N    <- length(V(gg)$name)
  meas <- matrix(0, nrow=N, ncol=3)

  for( i in 1:N ){
    sp <- as.numeric(igraph::shortest.paths(gg,i))
    sp <- sp[-i]
    sp <- sp[!sp == Inf]
    meas[i,1] <- min(sp)
    meas[i,2] <- round(mean(sp),3)
    meas[i,3] <- round(sd(sp),3)
  }

  return(meas)

}

#
# filterZeroDegree <- function( DIR, SUB ){
#
#   ft <- list()
#
#   for( d in 1:length(DIR) ){
#     tt <- read.table(sprintf("%s/random_%s_permuteDEGREE.csv",DIR[d],SUB),sep="\t",header=T)
#
#     Nc <- length(tt[1,])
#     tt <- tt[,2:Nc]
#
#     tmp <- apply(tt, 2, function(x) {ifelse(x == 0, NA, 1)})
#
#     ft[[d]] <- tmp
#
#     names(ft)[d] <- sprintf("RANDOM%d",d)
#
#   }
#
#   return(ft)
#
# }
#
#
# readRandomDataFiles <- function( DIR, SUB, CENT, C, FILTER , COLN ){
#
#   ranDF <- list()
#
#   for( d in 1:length(DIR) ){
#
#     if( C == 6 ){
#       tt <- read.table(sprintf("%s/random_MEAN_%s_permute%s.csv",DIR[d],SUB,CENT[C]),sep="\t",header=T)
#     } else {
#       tt <- read.table(sprintf("%s/random_%s_permute%s.csv",DIR[d],SUB,CENT[C]),sep="\t",header=T)
#     }
#
#
#     Nc <- length(tt[1,])
#
#     tt <- tt[,2:Nc]
#
#     if( !is.null(FILTER) ){
#       tt <- tt*FILTER[[d]]
#     }
#
#     ##test
#     tt  <- tt[,1]
#
#     tmp <- as.numeric(unlist(tt))
#
#     tmp <- tmp[!is.na(tmp)]
#
#     tmp <- data.frame(x=as.numeric(tmp),group=COLN[d])
#
#     ranDF[[d]] <- tmp
#
#     names(ranDF)[d] <- sprintf("RANDOM%d",d)
#
#     rm(tt,Nc,tmp)
#
#   }
#
#   return(ranDF)
#
# }
#
#
formatLogLogPlot <- function( X, GROUP ){

  X = as.vector(X)

  mm <- ecdf(X)

  df <- data.frame(x=sort(X),y=1-mm(sort(X)),group=GROUP)

  return(df)

}

filterLogLog <- function( df, xMAX, yMIN ){

  if( !is.null(xMAX) ){
    df <- df[ df$x <= xMAX, ]
  }

  if( !is.null(yMIN) ){
    df <- df[ df$y >= yMIN, ]
  }

  indx <- which(df$y==0)

  if( length(indx) != 0 ){
    df <- df[-indx,]
  }

  return(df)
}

##Calculate the Median absolute difference
MAD <- function( X ){

  X <- as.numeric(X)

  Xmd <- median(X)

  MAD <- median(abs(X-Xmd))

  return(MAD)
}

#' Calculate centrality measures for graph nodes.
#'
#' @param gg igraph object
#'
#' @return matrix with following columns:
#' * ID - vertex ID
#' * DEG - degree
#' * BET - betweenness
#' * CC - clustering coefficient
#' * SL - semilocal centrality
#' * mnSP - mean shortest path
#' * PR - page rank
#' * sdSP - standard deviation of the shortest path
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' m<-getCentralityMatrix(gg)
getCentralityMatrix<-function(gg){
  ID <- V(gg)$name
  N  <- length(ID)

  CN  <- c("ID","DEG","BET","CC","SL","mnSP","PR","sdSP")

  tmp <- matrix(0,nrow=N,ncol=length(CN))
  colnames(tmp) <- CN

  tmp[,1] <- ID
  tmp[,2] <- as.vector(igraph::degree(graph=gg))
  tmp[,3] <- as.character(round(betweenness(gg),3))
  tmp[,4] <- as.character(round(transitivity(gg,"local"),3))
  sl<- fSemilocal(gg)
  tmp[,5] <- as.character(round(sl,3))

  res <- calShorestPaths(gg)
  tmp[,6]  <- as.character(res[,2])
  tmp[,7]  <- as.character(round(as.vector(page.rank(graph=gg,vids=V(gg),directed=F,options=igraph.arpack.default)$vector),6))
  tmp[,8]  <- as.character(res[,3])

  return(tmp)

}

#' Add annotation to the vertex
#'
#' @param gg igraph object
#' @param m matrix of values to be applied as vertex attributes.
#'     matrix should contains column "ID" to map value to the vertex.
#'
#' @return modified igraph object
applpMatrixToGraph<-function(gg,m){
  ggm<-gg
  measures<-colnames(m)
  id.col<-which(measures=='ID')
  meas.col<-which(measures!='ID')
  for(i in meas.col){
    #remove previous annotation of that name
    #check does it really needed
    ggm<-removeVertexTerm(ggm,measures[i])
    idx<-match(V(gg)$name,m[,id.col])
    naid<-which(is.na(idx))
    if(length(naid)==0){
      ggm<-set.vertex.attribute(graph=ggm,
                                name=measures[i],
                                index = V(ggm),
                                value = m[idx,i])
    }else{
      gindex<-which(is.na(idx))
      ggm<-set.vertex.attribute(graph=ggm,
                                name=measures[i],
                                index = gindex,
                                value = m[idx[gindex],i])
    }
  }
  return(ggm)
}

#' Calculate centrality measures for graph nodes and save them as vertex
#' property.
#'
#' @param gg igraph object
#'
#' @return modified igraph object
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' ggm<-calcCentrality(gg)
calcCentrality<-function(gg){
  m<-getCentralityMatrix(gg)
  ggm<-applpMatrixToGraph(gg,m)
  return(ggm)
}

#get centrality measures for random graph
#' Title
#'
#' @param gg
#' @param type:
#' * gnp -- G(n,p) Erdos-Renyi model
#' * pa --  Barabasi-Albert model
#' * cgnp -- new random graph from a given graph by randomly adding/removing edges
#' * rw -- new random graph from a given graph by rewiring 25% of edges preserving the degree distribution
#'
#' @return
#' @export
#'
#' @examples
#' library(synaptome.db)
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
getRandomGraphCentrality<-function(gg,type=c('gnp','pa','cgnp','rw'),...){
  op<-options(warn= -1)
  type <- match.arg(type)
  nv<-vcount(gg)
  ne<-ecount(gg)
  prob<-(2*ne)/(nv*(nv-1))
    rg<-switch (type,
    gnp = getGNP(gg,...),
    pa  = getPA(gg,...),
    cgnp = sample_correlated_gnp(gg,corr=0.75,...),
    rw = rewire(gg,keeping_degseq(niter = 0.25*ne))
  )
  V(rg)$name<-V(gg)$name
  m<-getCentralityMatrix(rg)
  options(op)
  return(m)
}


getGNP<-function(gg,...){
  nv<-vcount(gg)
  ne<-ecount(gg)
  prob<-(2*ne)/(nv*(nv-1))
  g<-sample_gnp(nv,p=prob,...)
  return(g)
}

getPA<-function(gg,...){
  nv<-vcount(gg)
  pFit <- FitDegree( as.vector(igraph::degree(graph=gg)), Nsim=100, plot=FALSE )
  pwr <- pFit@alpha
  g<- sample_pa(nv,power=pwr,...)
  return(g)
}

#' Convert centrality matrix into
#'
#' @param m
#'
#' @return
#'
#' @examples
getGraphCentralityECDF<-function(m){
  idx<-which(colnames(m)!='ID')
  l<-list()
  for(i in 2:8){
    n<-colnames(m)[i]
    l[[n]]<-ecdf(as.numeric(m[,i]))
  }
  return(l)
}

#' Extracts particular measure from matrix and convert for distance calculation
#' by calcCentralityInternalDistances and calcCentralityExternalDistances
#' functions.
#'
#' @param m matrix of centrality measures as returned by getCentralityMatrix
#' @param nm name of the measure from m
#' @param keepOrder if FALSE valuess will be sorted
#'
#' @return
#'
#' @examples
getCM<-function(m,nm,keepOrder){
  v<-as.numeric(m[,which(colnames(m)==nm)])
  if(keepOrder){
    return(v)
  }else{
    return(sort(v,decreasing = FALSE,na.last=TRUE))
  }
}

#' Function calculates matrix of distances between elements of list
#'
#' @param l
#' @param keepOrder if FALSE valuess will be sorted
#' @param dist methods available from dist function
#'
#' @return
#'
#' @examples
calcCentralityInternalDistances<-function(l,keepOrder=FALSE,dist='euclidean'){
  CN  <- c("ID","DEG","BET","CC","SL","mnSP","PR","sdSP")
  resl<-list()
  for(i in 2:length(CN)){
    nm<-CN[i]
    res<-sapply(l,getCM,nm=nm,keepOrder=keepOrder)
    if(is.matrix(res)){
    resl[[nm]]<-as.vector(dist(t(res),method=dist))
    }
  }
  resm<-do.call(cbind,resl)
  return(resm)
}

#' Function calculates matrix of distances between elements of list and
#' the reference matrix
#'
#' @param m reference matrix
#' @param l list of permuted matrix
#' @param keepOrder if FALSE valuess will be sorted
#' @param dist methods available from dist function
#'
#' @return
#'
#' @examples
calcCentralityExternalDistances<-function(m,l,keepOrder=FALSE,dist='euclidean'){
  CN  <- c("ID","DEG","BET","CC","SL","mnSP","PR","sdSP")
  resl<-list()
  for(i in 2:length(CN)){
    nm<-CN[i]
    rm<-getCM(m,nm=nm,keepOrder=keepOrder)
    res<-sapply(l,getCM,nm=nm,keepOrder=keepOrder)
    if(is.matrix(res)){
    cmm<-cbind(rm,res)
    cmd<-as.matrix(dist(t(cmm),method=dist))
    resl[[nm]]<-as.vector(cmd[-1,1])
    }
  }
  resm<-do.call(cbind,resl)
  return(resm)
}

evalCentralitySignificance<-function(dmi,dme){
  nmi<-colnames(dmi)
  nme<-colnames(dme)
  nms<-intersect(nmi,nme)
  l<-list()
  for(nm in nms){
    mi<-dmi[,colnames(dmi)==nm]
    me<-dme[,colnames(dme)==nm]
    ks<-ks.test(mi,me)
    l[[nm]]<-list(ks=ks,
                  dt=data.frame(val=c(mi,me),
                                cl=factor(c(rep('perm',length(mi)),
                                     rep('graph',length(me))))))
  }
  return(l)
}



