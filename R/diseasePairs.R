qscore <- function(zz,FDR){

  LL <- FDR[FDR[,1] < as.numeric(zz),2]

  if( length(LL) != 0 ){ return(LL[end(LL)[1]]); }

  return(1)
}

##
# Calculate each diease-pair overlap/seperation on a selected
# synaptic PPI network models, based on analysis described in:
# Menche, J. et al. Uncovering disease-disease relationships through the incomplete interactome.
# Science, 347, (6224):1257601 (2015).
##

#source('../setUp.R')

#--- Find all Gene Disease Associations
#GDA <- V(gg)$TopOntoOVG

#Overlap of Disease A and B in the interactome
# GG   => igraph network
# GDA  => gda data for this graph
# disA => name of disease A
# disA => name of disease B
# OO   => minimum shorest paths for each gda, and each disease
diseaseOverlap <- function(GG, GDA, disA, disB, OO){

  #disease A genes
  IDS1  <- V(GG)$name[grepl(disA,GDA,fixed=T)]
  NIDS1 <- length(IDS1)

  #disease B genes
  IDS2  <- V(GG)$name[grepl(disB,GDA,fixed=T)]
  NIDS2 <- length(IDS2)

  #disease A given B
  paths  <- igraph::shortest.paths(GG,IDS1,IDS2,weights=NA)
  dsA    <- as.numeric(as.vector(apply(paths,1,min)))

  #disease B given A
  paths  <- igraph::shortest.paths(GG,IDS2,IDS1,weights=NA)
  dsB    <- as.numeric(as.vector(apply(paths,1,min)))

  #network-based separation between disease A and B
  dAB <- (sum(dsA)+sum(dsB))/(NIDS1+NIDS2)

  #network-based localisation of disease A
  indA <- which(colnames(OO)==disA)
  dA   <- mean(as.numeric(as.vector(OO[OO[,indA[1]]!=".",indA[1]])))

  #network-based localisation of disease B
  indB <- which(colnames(OO)==disB)
  dB   <- mean(as.numeric(as.vector(OO[OO[,indB[1]]!=".",indB[1]])))

  #overlap between disease A and B
  sAB = as.numeric(dAB) - (as.numeric(dA)+as.numeric(dB))/2

  return(sAB)

}

degree.binned.GDAs <- function(gg,GDA,dtype){

  deg  = degree(gg)
  bins = table(deg)
  map  = cbind(names(deg),as.vector(deg))
  map  = cbind(map,match(map[,2],names(bins)))

  nGDAs=length(dtype)

  for( i in 1:nGDAs){
    map=cbind(map,ifelse(grepl(dtype[i],GDA),1,0))
  }

  colnames(map) = c("EntrezID","Degree","Bin",dtype)

  return(map)

}

sample.deg.binned.GDA <- function(org.map,GDA){

  gda.indx = match(GDA,colnames(org.map))

  rnd.gene.set = NULL

  if( length(gda.indx) > 0 ){

    gda.set = as.vector(org.map[org.map[,gda.indx]==1,3])
    Nset    = length(gda.set)

    rnd.gene.set = rep(NA,Nset)
    map=org.map
    for( i in 1:Nset ){
      seq.map  = seq(1,dim(map)[1],1)
      rnd.indx = seq.map[!is.na(match(as.numeric(map[,3]),as.numeric(gda.set[i])))]
      if( length(rnd.indx) > 1 ){
        rnd.indx = as.numeric(sample(rnd.indx))[1]
      }
      if( length(rnd.indx) > 0 ){
        rnd.gene.set[i] = map[rnd.indx,1]
        map             = map[-rnd.indx,]
      }
    }
  }

  return(gdas=rnd.gene.set)

}


#' Get particular annotation from the graph and format it to the suitable
#' form.
#'
#' @param gg
#' @param name
#'
#' @return
#'
#' @examples
prepareGDA<-function(gg,name){
  gda<-get.vertex.attribute(gg,name)
  gda<-escapeAnnotation(gda)
  return(gda)
}

calcDiseasePairs<-function(gg,name,diseases=NULL,permute=c('none','random','binned')){
  permute<-match.arg(permute)
  gda<-prepareGDA(gg,name)
  if(is.null(diseases)){
    diseases<-getAnnotationList(gda,sort='freq')
  }else{
    remove<-c()
    diseases<-escapeAnnotation(diseases)
    for(d in 1:length(diseases)){
      if(!any(grepl(diseases[d],gda))){
        remove<-c(remove,d)
      }
    }
    if(length(remove)>0){
      diseases<-diseases[-remove]
    }
  }
  if(permute=='binned'){
    map <- degree.binned.GDAs(gg,gda,diseases)
  }
  res           <- matrix(0 ,ncol=4, nrow=length(diseases))
  colnames(res) <- c("Disease","N","mean_ds","SD_ds")
  res[,1]       <- diseases


  #--- store minimum shorest paths for each gda, and each disease
  oo <- matrix(".",nrow=NN,ncol=(length(diseases)+2))
  colnames(oo) <- c("Gene.ID","Gene.Name",diseases)
  oo[,1]       <- V(gg)$name[gda !=""]
  oo[,2]       <- V(gg)$GeneName[gda !=""]

  ##--- loop over each disease
  for( d in 1:length(diseases) ){

    IDS <- V(gg)$name[grepl(diseases[d],gda,fixed=T)]
    N   <- length(IDS)

    if(permute=='random'){
      ## permute the N GDA's relative to all gene ids
      IDS <- permute(GNS, NN) #case
    }else if(permute=='binned'){
      IDS <- sample.deg.binned.GDA(map,diseases[d])
    }

    ## for each gda, find the minimum shortest path to next gda (of the same disease)
    XX=igraph::shortest.paths(gg,IDS,IDS,weights=NA)
    diag(XX)       = NA
    ds             = apply(XX,1,min,na.rm=T)
    indX           = match(names(ds),oo[,1])
    oo[indX,(2+d)] = as.vector(ds)

    res[d,2]       = as.character(N)
    res[d,3]       = as.character(mean(ds))
    res[d,4]       = as.character(sd(ds))

  }

  DAB <- matrix(".",ncol=length(diseases),nrow=length(diseases))
  colnames(DAB) <- diseases
  rownames(DAB) <- diseases

  #--- NOTE ---#
  # DAB is bound by -dmax <= DAB <= dmax
  # where dmax denotes the diameter of the network
  # dmax <- diameter(gg,directed=F)
  #------------#

  ##--- calculate disease-disease overlap
  for( i in 1:length(diseases) ){
    for( j in i:length(diseases) ){

      DAB[i,j] <- 0

      if( i != j ){
        DAB[i,j] <- diseaseOverlap(gg,gda,rownames(DAB)[i],colnames(DAB)[j],oo)
      }

    }
  }

  return(list(disease_separation=DAB,gene_disease_separation=oo,disease_localisation=res))
}

runPermDisease<-function(gg,name,diseases=NULL,Nperm=100){
  res<-lapply(1:Nperm,calcDiseasePairs,gg=gg,name=name,diseases=diseases,perm=TRUE)
}
