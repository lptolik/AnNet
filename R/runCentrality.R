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
  tmp[,5] <- as.character(round(Semilocal(gg),3))

  res <- as.matrix(calShorestPaths(gg))
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
#'
#' @return
#' @export
#'
#' @examples
getRandomGraphCentrality<-function(gg,type=c('gnp','pa'),...){
  nv<-vcount(gg)
  ne<-ecount(gg)
  rg<-switch (type,
    gnp = sample_gnp(nv,...),
    pa  = sample_pa(nv,...)
  )
  m<-getCentralityMatrix(rg)
  return(m)
}
