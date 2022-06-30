#---Get all edges internal to a community
intraEdges <- function(GG, ALG, CC, INTRA=NULL, INTER=NULL){

  intra = NULL #edges in the community CC
  inter = NULL #edges going out from community CC

  if( !is.null(igraph::get.vertex.attribute(GG,ALG)) ){

    coms <- get.vertex.attribute(GG,ALG)

    if( length(which(coms == CC)) != 0 ){

      ed_cc = E(GG)[inc(coms == CC)]

      all_edges_m <- get.edges(GG, ed_cc) #matrix representation

      inter = (ed_cc[!(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])

      intra = (ed_cc[(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])

    }
  }

  if( INTRA==TRUE && !is.null(intra) ){
    intra_m = get.edges(GG,intra)
    intra   = cbind(V(GG)$name[intra_m[,1]],V(GG)$name[intra_m[,2]])
    return(intra)
  }

  if( INTER==TRUE && !is.null(inter) ){
    inter_m = get.edges(GG,inter)
    inter   = cbind(V(GG)$name[inter_m[,1]],V(GG)$name[inter_m[,2]])
    return(inter)
  }

  return(NULL)

}

intraEdgesM <- function(GG, mem, CC, INTRA=NULL, INTER=NULL){

  intra = NULL #edges in the community CC
  inter = NULL #edges going out from community CC

  idx<- (mem$membership == CC)
    if( length(which(idx)) != 0 ){

      ed_cc = E(GG)[inc(idx)]

      all_edges_m <- get.edges(GG, ed_cc) #matrix representation

      inter = (ed_cc[!(all_edges_m[, 1] %in% V(GG)[idx] & all_edges_m[, 2] %in% V(GG)[idx])])

      intra = (ed_cc[(all_edges_m[, 1] %in% V(GG)[idx] & all_edges_m[, 2] %in% V(GG)[idx])])

    }

  if( INTRA==TRUE && !is.null(intra) && length(intra)>0 ){
    intra_m = get.edges(GG,intra)
    intra   = cbind(V(GG)$name[intra_m[,1]],V(GG)$name[intra_m[,2]])
    return(intra)
  }

  if( INTER==TRUE && !is.null(inter) && length(inter)>0 ){
    inter_m = get.edges(GG,inter)
    inter   = cbind(V(GG)$name[inter_m[,1]],V(GG)$name[inter_m[,2]])
    return(inter)
  }

  return(NULL)

}

#' Return induced subgraph for cluster
#'
#' Function takes graph \code{gg}, membership data.frame \code{mem} and
#' ID of the cluster in it, creates induced subgraph and returned it.
#'
#' @param clID
#' @param gg
#' @param mem
#'
#' @return induced subgraph as igraph object
#' @export
#'
#' @examples
getClusterSubgraphByID<-function(clID,gg,mem){
  idx<-which(mem$membership==clID)
  sg<-induced_subgraph(gg,V(gg)[idx],impl = "auto")
  return(sg)
}

#' Calculate layout based upon membership
#'
#' @param gg
#' @param mem
#'
#' @return
#' @export
#'
#' @examples
layoutByCluster<-function(gg,mem,layout=layout_with_kk){
  Cn<-table(mem$membership)
  sgraphs<-lapply(names(Cn),getClusterSubgraphByID,gg=gg,mem=mem)
  layouts <- lapply(sgraphs, layout)
  lay <- merge_coords(sgraphs, layouts)
  ug <- disjoint_union(sgraphs)
  idx<-match(V(gg)$name,V(ug)$name)
  lay<-lay[idx,]
  return(lay)
}

#' Calculate two-level layout from recluster matrix
#'
#' Takes results of recluster and apply \code{layoutByCluster} to each
#'
#' @param gg
#' @param remem
#' @param layout
#'
#' @return
#' @export
#'
#' @examples
layoutByRecluster<-function(gg,remem,layout=layout_with_kk){
  Cn<-table(remem$membership)
  glist<-list()
  laylist<-list()
  for(i in 1: length(Cn)){
    sg<-getClusterSubgraphByID(names(Cn)[i],gg,remem$membership)
    mem1<-remem[remem$membership==names(Cn)[i],c('names','recluster')]
    names(mem1)<-c('names','membership')
    if(length(table(mem1$membership))>1){
      lay<-layoutByCluster(sg,mem1,layout)
    }else{
      lay<-layout(sg)
    }
    glist[[i]]<-sg
    laylist[[i]]<-lay
  }
  ug <- disjoint_union(glist)
  lay <- merge_coords(glist, laylist)
  idx<-match(V(gg)$name,V(ug)$name)
  layF<-lay[idx,]
  return(layF)
}

#' Recluster graph
#'
#' Function takes graph \code{gg} and its membership matrix \code{mem}
#' as returned \code{calcMembership} and apply clustering algorithm \code{alg}
#' to all clusters larger than \code{CnMAX}
#'
#' @param gg
#' @param mem
#' @param alg
#' @param CnMAX
#'
#' @return
#' @export
#'
#' @examples
calcReclusterMatrix<-function(gg,mem,alg,CnMAX,keepSplit=FALSE){

  if(is.matrix(mem)){
    mem<-as.data.frame(mem)
  }

  if(!all(c('names','membership')%in%names(mem))){
    stop("mem suppose to have columns 'names' and 'membership'")
  }
    #--- algorithm clustering 1
    ALG1 <- mem

    Cn <- table(mem$membership)
    Cnc <- Cn[Cn > CnMAX]
    cc <- names(Cnc)


    RES <- list()
    k=1
    for( i in 1:length(cc) ){

      edCC = intraEdgesM(gg, mem, cc[i], INTRA=TRUE)

      if( !is.null(edCC) ){

        ggLCC    <- graph_from_data_frame(d=edCC, directed=F)
        res <- getClustering(ggLCC,alg)
        oo       <- data.frame(names=res$names, membership=res$membership)
        if(dim(oo)[1]< Cnc[i]){
          cmem<-mem[mem$membership==cc[i]]
          singidx<-which(!cmem$names %in% oo$names)
          singletones <- data.frame(names=cmem$names[singidx],
                                    membership=max(oo$membership)+
                                      1:length(singidx))
          oo<-rbind(oo,singletones)
        }

        RES[[k]]      <- oo
        names(RES)[k] <- cc[i]
        k=k+1
      }

    }#for


    if( length(RES) == 0 ){ return(NULL) }


    #--- algorithm clustering 2
    ALG2     <- mem
    ALG2$split <- rep(-1, dim(ALG1)[1])
    indx     <- match(ALG2$membership,cc)
    indx     <- ifelse(is.na(indx),TRUE, FALSE)
    ALG2$split <- ifelse(indx, ALG2$membership, ALG2$split)

    CCmax = max(as.numeric(ALG2$split))

    for( i in 1:length(cc) ){

      temp     <- RES[[i]]
      temp$membership <- temp$membership + CCmax

      indx <- match(ALG2$names,temp$names)

      ALG2$split <- ifelse(is.na(indx),ALG2$split,temp$membership[indx])

      CCmax = max(as.numeric(ALG2$split))

    }

    #---reorder ALG2$split
    N = length(V(gg));

    temp    <- rep(-1, N)
    counter <- min(as.numeric(ALG2$split))
    Knew    <- 1;
    Kmax    <- max(as.numeric(ALG2$split))

    while( counter <= Kmax ){

      found=FALSE;

      for(v in 1:N ){
        if( as.numeric(ALG2$split[v]) == counter ){
          temp[v] = Knew;
          found=TRUE;
        }
      }

      if(found) Knew=Knew+1;

      counter=counter+1;
    }

    #---final
    ALG3 <- cbind(ALG2, data.frame(recluster=temp))
    if(!keepSplit){
      ALG3<-ALG3[,grep('split',names(ALG3),invert = TRUE)]
    }
    return(ALG3)


}

recluster <- function( GG, ALGN, CnMAX ){

  if( !is.null(igraph::get.vertex.attribute(GG,ALGN)) ){


    #--- algorithm clustering 1
    ALG1 <- get.vertex.attribute(GG,ALGN,V(GG))
    ALG1 <- cbind(V(GG)$name, ALG1)

    Cn <- table(as.numeric(ALG1[,2]))
    cc <- names(Cn)[Cn > CnMAX]


    RES <- list()
    k=1
    for( i in 1:length(cc) ){

      edCC = intraEdges(GG, ALG, cc[i], INTRA=TRUE)

      if( !is.null(edCC) ){

        ggLCC    <- graph_from_data_frame(d=edCC, directed=F)
        res <- getClustering(ggLCC,alg)
        oo       <- cbind(res$names, res$membership)

        RES[[k]]      <- oo
        names(RES)[k] <- cc[i]
        k=k+1
      }

    }#for


    if( length(RES) == 0 ){ return(NULL) }


    #--- algorithm clustering 2
    ALG2     <- cbind(ALG1, rep(-1, length(ALG1[,1])))
    indx     <- match(ALG2[,2],cc)
    indx     <- ifelse(is.na(indx),TRUE, FALSE)
    ALG2[,3] <- ifelse(indx, ALG2[,2], ALG2[,3])

    CCmax = max(as.numeric(ALG2[,3]))

    for( i in 1:length(cc) ){

      temp     <- RES[[i]]
      temp[,2] <- as.numeric(temp[,2]) + CCmax

      indx <- match(ALG2[,1],temp[,1])
      indx <- temp[indx,2]

      ALG2[,3] = ifelse(is.na(indx),ALG2[,3],indx)

      CCmax = max(as.numeric(ALG2[,3]))

    }

    #---reorder ALG2[,3]
    N = length(V(GG));

    temp    <- rep(-1, N)
    counter <- min(as.numeric(ALG2[,3]))
    Knew    <- 1;
    Kmax    <- max(as.numeric(ALG2[,3]))

    while( counter <= Kmax ){

      found=FALSE;

      for(v in 1:N ){
        if( as.numeric(ALG2[v,3]) == counter ){
          temp[v] = Knew;
          found=TRUE;
        }
      }

      if(found) Knew=Knew+1;

      counter=counter+1;
    }

    #---final
    ALG3 <- cbind(ALG2, temp)
    return(ALG3)
  }

  return(NULL)

}

#aa    <- as.numeric(args[1]) #
#type  <- as.numeric(args[2]) #edges=>1 or nodes=>2  to mask
#mask  <- as.numeric(args[3]) #number of edges/nodes to mask
#Cnmin <- as.numeric(args[4]) #Cn min for Spectral algorithm
#Cnmax <- as.numeric(args[5]) #Cn max for reclustering algorithms

#' Perturbe graph and calculate its clustering
#'
#' @param gg graph
#' @param mask percentage of elements to perturbe
#' @param alg clustering alg.
#' @param type edges=>1 or nodes=>2  to mask
#' @param reclust logical to decide wether to invoke reclustering
#' @param Cnmin Cn min for Spectral algorithm
#' @param Cnmax Cn max for reclustering algorithms
#'
#' @return
#' @export
#'
#' @examples
sampleGraphClust<-function(gg,mask,alg,type,reclust=FALSE,Cnmin=-1,Cnmax=10){
  IDS <- V(gg)$name;
  ids <- V(gg)$name;

  #---subsampling scheme
  if( type == 1 ){

    nr  <- ceiling( length(E(gg))*(mask/100) )
    ggM <- delete_edges(gg,sample(E(gg),nr))

  }

  if( type == 2 ){

    nr  <- ceiling( length(V(gg))*(mask/100) )
    ggM <- delete_vertices(gg,sample(V(gg),nr))

  }


  #---Find Largest CC
  ggLCC <- findLCC(ggM)
  #---


  #---build consensus file
  cc       <- matrix(-1, ncol=3, nrow=length(V(gg)))
  cc[,1]   <- V(gg)$name
  cc[,2]   <- ifelse(cc[,1] %in% V(ggLCC)$name,cc[,1],-1)


  if( Cnmin > 0 ){
    Cnmin = floor( (Cnmin*length(V(ggLCC)))/100 )
  } else {
    Cnmin = 1;
  }
  if( Cnmax > 0 ){
    Cnmax = floor( (Cnmax*length(V(ggLCC)))/100 )
  } else {
    #---default is 10% of network size
    Cnmax = floor( (10*length(V(ggLCC)))/100 )
  }
  cl<-getClustering(ggLCC,alg)
  if(reclust){
    ggLCC = igraph::set.vertex.attribute(ggLCC,alg,V(ggLCC),cl$membership)

    oo = recluster( ggLCC, alg, Cnmax )

    if( !is.null(oo) ){
      cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }


  }else{
  cc[,3]   <- ifelse(cc[,2] %in% cl$names,cl$membership,-1)
  }
  return(cc)
}
