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

getEntropy<-function(gg){
  #--- initial entropy rate
  V    <- length(V(gg))
  E    <- length(E(gg))
  ki   <- as.vector(igraph::degree(gg))
  Kbar <- mean(ki)

  #--- get adjacency matrix for graph
  A    <- get.adjacency(gg)

  #--- get leading eigenvalue and vector
  R     <- eigen(A)
  Rindx <- which.max(R$values)
  gamma <- R$values[Rindx]
  nu    <- R$vectors[,Rindx]


  #--- calculate max entropy rate, maxSr
  Pij   <- (A * nu) / (gamma * nu)
  Pi    <- ki/(2*E)
  maxSr <- sum(as.vector(Pi * apply(Pij,1,maxLSi,BASE=0)))


  #--- calculate initial configuration
  Norm <- as.numeric(V*Kbar) #as.numeric(2*E)
  SRo  <- as.numeric(1/Norm)*sum(ki*log(ki))


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
  return(SRprime)
}
