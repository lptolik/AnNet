## build consensus matrix "C" from a set
## of randomised results files from a
## clustering algorithm
#' Build consensus matrix from set of files.
#'
#' @param Dir
#' @param file.name
#' @param skip
#' @param sep
#'
#' @return
buildConsensusMatFromFiles <- function(Dir,file.name,skip=1,sep="\t"){

  subdirs  = list.files(path=Dir,pattern=file.name);
  nstudies = length(subdirs);

  N        = NULL
  I        = NULL
  M        = NULL
  C        = NULL
  initM    = TRUE

  NJobs    = 0
  max_com  = 0
  min_com  = 500
  for( s in 1:nstudies ){
    filein = sprintf("%s/%s/%s",Dir,subdirs[s]);
    if( file.exists(filein) && file.info(filein)$size!=0 ){
      tb = read.delim(filein,
                      skip=skip,
                      header=F,
                      sep=sep);
      ## make sure node id == -1 if node com == -1
      indx = tb[,3] == -1
      tb[indx,2] = -1
      if(initM){
        N     = dim(tb)[1]
        temp  = matrix(0,nrow=N,ncol=N)
        I     = temp
        M     = temp
        initM = FALSE
        rm(temp)
      }
      k.coms    = tb[tb[,3] != -1,3]
      k.max     = max(k.coms,na.rm=T)
      k.min     = min(k.coms,na.rm=T)
      if( k.max > max_com   ){ max_com   = k.max; }
      if( k.min < min_com   ){ min_com   = k.min; }
      study = calculateConsensusMat( data=tb )
      I = I + study$I
      M = M + study$M
      NJobs = NJobs + 1;
      rm(tb,study)
    }
  }

  ## the consensus matrix
  if( !is.null(N) ){
    C = do.call( cbind, lapply(1:N, function(s) matrixDiv(M[,s],I[,s])))
  }
  return(C)
}

#' Build consensus matrix from the list of clusterings
#'
#' @param lcc list of clustering matrices obtained from the
#' \code{\link{sampleGraph}}
#'
#' @return consensus matrix
buildConsensusMatrix<-function(lcc){
  N        = NULL
  I        = NULL
  M        = NULL
  C        = NULL
  initM    = TRUE

  NJobs    = 0
  max_com  = 0
  min_com  = 500
  for(i in 1:length(lcc)){
    tb<-lcc[[i]]
    ## make sure node id == -1 if node com == -1
    indx = tb[,3] == -1
    tb[indx,2] = -1

    if(initM){
      N     = dim(tb)[1]
      temp  = matrix(0,nrow=N,ncol=N)
      I     = temp
      M     = temp
      initM = FALSE
      rm(temp)
    }
    study = calculateConsensusMat( data=tb )
    I = I + study$I
    M = M + study$M
    NJobs = NJobs + 1;
    rm(tb,study)
  }
  if( !is.null(N) ){
    C = do.call( cbind, lapply(1:N, function(s) matrixDiv(M[,s],I[,s])))
  }
  return(C)
}

## divide columns of matrix X by Y
matrixDiv <- function(x,y){

  N    = length(x)
  res  = rep(0,N)
  indx = y != 0
  if( sum(indx) != 0 ){
    res[indx] = x[indx]/y[indx]
  }

  return(res)
}

#' Function to build consensus matrix in memory
#'
#' @param gg graph to perturb
#' @param N number of perturbation steps
#' @param mask percentage of elements to perturbe
#' @param alg clustering alg.
#' @param type edges=>1 or nodes=>2  to mask
#' @param reclust logical to decide wether to invoke reclustering
#' @param Cnmin Cn min for Spectral algorithm
#' @param Cnmax Cn max for reclustering algorithms
#'
#' @return consensus matrix of Nvert X Nvert
#' @export
#'
#' @examples
#' data(karate,package='igraphdata')
#' alg<-'louvain'
#' gg<-calcClustering(karate,alg = alg)
#' conmat<-makeConsensusMatrix(gg,N=100,mask = 10,alg = alg,type = 2)
#' dim(conmat)
makeConsensusMatrix<-function(gg,N=500,mask=20,alg,type,
                              reclust=FALSE,Cnmin=-1,Cnmax=10){
  lcc<-lapply(1:N, function(.x)sampleGraphClust(gg=gg,
                                            mask=mask,
                                            alg=alg,
                                            type=type,
                                            reclust=reclust,
                                            Cnmin=Cnmin,
                                            Cnmax=Cnmax))
  mm<-buildConsensusMatrix(lcc)
  colnames(mm)<-V(gg)$name
  rownames(mm)<-V(gg)$name
  return(mm)
}

## calculate the identity matrix "I" and tally matrix "M"
## given each node pairs community assignment
calculateConsensusMat <- function( data=NULL ){

  I = NULL
  M = NULL

  if( !is.null(data) ){

    N     = dim(data)[1]
    tempI = matrix(0,nrow=N,ncol=N)
    tempM = matrix(0,nrow=N,ncol=N)

    for( i in 1:N ){
      comi = as.numeric(data[i,3])
      keyi = data[i,2]
      jj   = seq(i,N,1)
      if( keyi != '-1' ){

        comj = as.numeric(data[jj,3])
        keyj = data[jj,2]

        ## I
        indxJ = jj[keyj!='-1']
        Nindx = length(indxJ)

        tempI[i,indxJ] = as.numeric(rep(1,Nindx))
        tempI[i,i]     = as.numeric(0.5)

        ## M
        indxC = jj[comj != -1 & comi == comj]
        Nindx = length(indxC)

        tempM[i,indxC] = as.numeric(rep(1,Nindx))
        tempM[i,i]     = as.numeric(0.5)

      }
    }

    M = tempM + t(tempM)
    I = tempI + t(tempI)

    rm(tempM,tempI)

  }

  return(list(M=M,I=I))

}



