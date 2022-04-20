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


recluster <- function( GG, ALGN, CnMAX, CnMIN=1 ){

  if( !is.null(igraph::get.vertex.attribute(GG,ALGN)) ){


    #--- algorithm clustering 1
    ALG1 <- get.vertex.attribute(GG,ALGN,V(GG))
    ALG1 <- cbind(V(GG)$name, ALG1)

    Cn <- table(as.numeric(ALG1[,2]))
    cc <- names(Cn)[Cn > CnMAX]


    RES <- list()
    k=1
    for( i in 1:length(cc) ){

      edCC = intraEdges(GG, ALGN, cc[i], INTRA=TRUE)
      oo       <- cbind(res$names, res$membership)

      if( !is.null(edCC) ){

        ggLCC    <- graph_from_data_frame(d=edCC, directed=F)
        res <- getClustering(ggLCC,alg)

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
    ggLCC = igraph::set.vertex.attribute(ggLCC,alg,V(ggLCC),louvain$membership)

    oo = recluster( ggLCC, alg, Cnmax )

    if( !is.null(oo) ){
      cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }


  }else{
  cc[,3]   <- ifelse(cc[,2] %in% cl$names,cl$membership,-1)
  }
  return(cc)
}
