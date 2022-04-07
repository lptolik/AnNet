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
#' @export
#'
#' @examples
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

  cat("build consensus matrix...")

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
  cat("done!\n")
  cat("Number of good Jobs: ", NJobs, "\n")
  cat("max.com            : ", max_com, "\n")
  cat("min.com            : ", min_com, "\n")
  ##cat("avg.com            : ", max_com/NJobs, "\n")

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
      keyi = as.numeric(data[i,2])
      jj   = seq(i,N,1)
      if( keyi != -1 ){

        comj = as.numeric(data[jj,3])
        keyj = as.numeric(data[jj,2])

        ## I
        indxJ = jj[keyj!=-1]
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



