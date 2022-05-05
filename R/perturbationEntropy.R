maxLSi <- function( XX, BASE=0 ){

  XX  <- as.vector(as.numeric(XX))

  XXo <- XX[XX != 0]
  if( BASE == 2 ){
    return ( -sum( XXo * log2(XXo)) )
  }

  if( BASE == 10 ){
    return ( -sum( XXo * log10(XXo)) )
  }

  if( BASE == 0 ){
    return ( -sum( XXo * log(XXo)) )
  }

}

#' Calculate parameters for Entropy plot and calculations.
#'
#' @param gg igroph object
#'
#' @return list with values of maxSr and SRo
#' @export
#' @import RSpectra
#'
#' @examples
getEntropyRate<-function(gg){
  V    <- length(V(gg))
  E    <- length(E(gg))
  ki   <- as.vector(igraph::degree(gg))
  Kbar <- mean(ki)

  #--- get adjacency matrix for graph
  A    <- get.adjacency(gg)

  #--- get leading eigenvalue and vector
  #R     <- eigen(A)
  #Rindx <- which.max(R$values)
  #gamma <- R$values[Rindx]
  #nu    <- R$vectors[,Rindx]
  R<-RSpectra::eigs(A,1)
  gamma <- R$values[1]
  nu    <- R$vectors[,1]


  #--- calculate max entropy rate, maxSr
  Pij   <- (A * nu) / (gamma * nu)
  Pi    <- ki/(2*E)
  maxSr <- sum(as.vector(Pi * apply(Pij,1,maxLSi,BASE=0)))


  #--- calculate initial configuration
  Norm <- as.numeric(V*Kbar) #as.numeric(2*E)
  SRo  <- as.numeric(1/Norm)*sum(ki*log(ki))

  return(list(maxSr=maxSr,SRo=SRo))
}

#' Calculates perturbation entropy
#'
#' @param gg igraph object
#'
#' @return
#' @export
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' gg<-annotateGeneNames(gg)
#' e<- getEntropy(gg)
getEntropy<-function(gg,maxSr=NULL){
  if(!"GeneName"%in%vertex_attr_names(gg)){
    V(gg)$GeneName<-V(gg)$name
  }
  #--- initial entropy rate
  V    <- length(V(gg))
  E    <- length(E(gg))
  ki   <- as.vector(igraph::degree(gg))
  Kbar <- mean(ki)

  #--- get adjacency matrix for graph
  A    <- get.adjacency(gg)

  if(is.null(maxSr)){
    par<-getEntropyRate(gg)
    maxSr<-par$maxSr
  }

  #--- perturbated each PPI node/gene

  #--- expression values
  xx    <- vector(length=2)
  xx[1] <- 2  #active
  xx[2] <- 16 #inactive

  #--- perturbation expression values
  lambda    <- vector(length=2)
  lambda[1] <- 14   #active
  lambda[2] <- -14  #inactive

  #--- Norm for PI'
  NORM      <- vector(length=2)
  NORM[1]   <- 0
  NORM[2]   <- 0

  SRprime <- cbind(V(gg)$name, V(gg)$GeneName, ki, rep("",V), rep("",V))

  for( v in 1:V ){

    #--- name of gene to perturb
    GN     <- as.character(SRprime[v,2])
    GNindx <- which(V(gg)$GeneName==GN)

    #--- PI'
    PIprime <- cbind( rep("",V), rep("",V) )

    #--- LS'
    LSprime <- cbind( rep("",V), rep("",V) )

    #--- reset NORM
    NORM[1] = 0; NORM[2] = 0;

    #--- calculate norm for PI'
    for( s in 1:length(lambda) ){
      X               <- rep(xx[s], V)
      X[GNindx[1]]    <- X[GNindx[1]] + lambda[s]
      NORM[s]         <- X %*% A %*% X
    }


    #--- find all neighors to v, i.e. N(v)
    Nv <- V(gg)$name[neighbors(gg,GNindx,mode="all")]

    oo <- cbind( ki, !(V(gg)$name %in% Nv) )


    #--- PI' when v is not N(v)
    for( s in 1:length(lambda) ){
      PIprime[,s] <- ifelse(oo[,2] == 1, (1/NORM[s] * xx[s] * xx[s] * as.numeric(oo[,1])), ".")
    }

    #--- PI' when v is v
    for( s in 1:length(lambda) ){

      X   <- as.numeric(xx[s])
      lam <- as.numeric(lambda[s])
      DEG <- as.numeric(oo[GNindx[1],1])

      PIprime[GNindx[1],s] <- ((X + lam) * DEG * X) / NORM[s]

    }

    #--- PI' when v is N(v)
    for( s in 1:length(lambda) ){
      PIprime[,s] <- ifelse(oo[,2] == 0, (1/NORM[s] * xx[s] * ( xx[s] + lambda[s] + (as.numeric(oo[,1]) - 1) * xx[s])),PIprime[,s])
    }


    #--- LS' when v is not N(v)
    for( s in 1:length(lambda) ){
      X <- as.numeric(xx[s])
      LSprime[,s] <- ifelse(oo[,2] == 1, (-log(X) + log(X*as.numeric(oo[,1]))),".")
    }

    #--- LS' when v is N(v)
    Ni <- grep(0,oo[,2])

    for( i in 1:length(Ni) ){

      DEGi <- as.numeric(oo[Ni[i],1])
      SUM  <- DEGi-1

      for( s in 1:length(lambda) ){

        X   <- as.numeric(xx[s])
        lam <- as.numeric(lambda[s])

        dem <- X + lam + (DEGi -1) * X

        pij <- X / dem
        pi1 <- (X + lam) / dem

        #cat("v",v, ": i ", i, ": pij ", pij, ": pi1 ", pi1,"\n")

        LSi <- pij * log(pij)
        LSi <- - SUM * LSi - pi1 * log(pi1)

        LSprime[Ni[i],s]  <- as.character(LSi)

      }

    }

    SRprime[v,4] <- sum( as.numeric(PIprime[,1]) * as.numeric(LSprime[,1]) )
    SRprime[v,5] <- sum( as.numeric(PIprime[,2]) * as.numeric(LSprime[,2]) )

  }
  # colnames(SRprime) <- c("ENTREZ.ID","GENE.NAME","DEGREE","UP","DOWN")
  # SRprime <- as.data.frame(SRprime)
  # SRprime$DEGREE<-as.numeric(SRprime$DEGREE)
  # SRprime$UP <- as.numeric(SRprime$UP)/maxSr
  # SRprime$DOWN <- as.numeric(SRprime$DOWN)/maxSr

  SRprime[,4] <- as.numeric(SRprime[,4])/maxSr
  SRprime[,5] <- as.numeric(SRprime[,5])/maxSr

  colnames(SRprime) <- c("ENTREZ.ID","GENE.NAME","DEGREE","UP","DOWN")


  return(SRprime)
}

getEntropyOverExpressed<-function(SRprime,perc=1){
  #--- Bottom 1% UP, i.e. OVER-EXPRESSED
  V <- dim(SRprime)[1]
  XX  <- as.numeric(SRprime[,4])
  MIN <- min(XX)
  MAX <- max(XX)
  XX2 <- (XX-MIN)/(MAX-MIN)
  oo2 <- cbind(SRprime[,2],XX2)
  oo2 <- oo2[order(as.numeric(oo2[,2])),]
  ii  <- floor(perc/100 * V)
  GN  <- oo2[1:ii,1]
  DF3 <- SRprime[match(GN,SRprime[,2]),c(1,2,3,4)]
  DF3 <- cbind(DF3,rep(paste0(perc,"%"),length(GN)))
  return(DF3)
}

#' Plot entropy values
#'
#' @param SRprime
#' @param subTIT
#' @param SRo
#' @param maxSr
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
plotEntropy<-function(SRprime,subTIT='Entropy',SRo=NULL,maxSr=NULL){
  colours <- c('lawngreen','firebrick2')

  V <- dim(SRprime)[1]
  DF1 <- SRprime[,c(1,2,3,4)]
  DF1 <- cbind(DF1,rep("SR_UP",length(SRprime[,1])))

  DF2 <- SRprime[,c(1,2,3,5)]
  DF2 <- cbind(DF2,rep("SR_DOWN",length(SRprime[,1])))

  DF  <- rbind(DF1,DF2)
  colnames(DF) <- c("ENTREZ.ID","GENE.NAME","DEGREE","SR","GROUP")

  DF <- as.data.frame(DF)
  DF$DEGREE<-as.numeric(DF$DEGREE)
  DF$SR<-as.numeric(DF$SR)
  DF$GROUP <- as.factor(DF$GROUP)

  gplot <- ggplot(DF,aes(x=log(DEGREE),y=SR, colour=GROUP) )+
    geom_point()+
    labs(x="log(k)",y="SR",title=subTIT)+
    guides(color=guide_legend(override.aes=list(fill=NA,size=4)),
           fill  = FALSE,
           group = FALSE,
           alpha = FALSE)+
    theme(
      axis.title.x=element_text(face="bold",size=rel(2)),
      axis.text.x =element_text(face="bold",size=rel(2)),
      axis.title.y=element_text(face="bold",size=rel(2)),
      axis.text.y =element_text(face="bold",size=rel(2)),
      legend.title=element_text(face="bold",size=rel(1.5)),
      legend.text=element_text(face="bold",size=rel(1.5)),
      legend.position="top",
      legend.key=element_blank())+
    theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
          panel.grid.minor = element_line(colour="grey40",size=0.1),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype="solid",fill=NA))+
    scale_color_manual("",breaks=levels(DF$GROUP),values=c(colours))#+
  if(!is.null(SRo) && !is.null(maxSr)){
    gplot <- gplot + geom_hline(yintercept=SRo/maxSr,colour="black",size=2,linetype=2,show.legend=F)
  }
  #geom_hline(yintercept=SRo/maxSr,colour="grey40",size=2,linetype=2,show.legend=F)
  return(gplot)
}
