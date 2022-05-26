#' @export
qscore <- function(zz,FDR){

  LL <- FDR[FDR[,1] < as.numeric(zz),2]

  if( length(LL) != 0 ){ return(LL[end(LL)[1]]); }

  return(1)
}

#' Title
#'
#' @param GNS
#' @param N
#'
#' @return
#' @export
#'
#' @examples
permute <- function(GNS, N){

  temp <- sample(GNS,N,replace=F)

  return(temp)

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
#' @export
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

#' @export
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

#' @export
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
#' @export
#'
#' @examples
prepareGDA<-function(gg,name){
  gda<-get.vertex.attribute(gg,name)
  gda<-escapeAnnotation(gda)
  return(gda)
}

#' Title
#'
#' @param gg
#' @param name
#' @param diseases
#' @param permute
#'
#' @return
#' @export
#'
#' @examples
calcDiseasePairs<-function(gg,name,diseases=NULL,permute=c('none','random','binned')){
  permute<-match.arg(permute)
  gda<-prepareGDA(gg,name)
  NN  <- length(which(gda!=""))
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
      IDS <- AnNet::permute( oo[,1], N) #case
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

  DAB <- matrix(NA,ncol=length(diseases),nrow=length(diseases))
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

#' Title
#'
#' @param gg
#' @param name
#' @param diseases
#' @param Nperm
#'
#' @return
#' @export
#'
#' @examples
runPermDisease<-function(gg,name,diseases=NULL,Nperm=100){
  resD<-calcDiseasePairs(gg=gg,name=name,diseases=diseases,permute = 'random')
  ds<-resD$gene_disease_separation
  loc<-resD$disease_localisation
  resL<-lapply(1:Nperm,function(.x)calcDiseasePairs(gg=gg,name=name,diseases=diseases,permute='random'))
  resA<-sapply(resL,function(.x).x$disease_separation,simplify = "array")
  resM<-apply(resA,c(1,2),mean)
  resS<-apply(resA,c(1,2),sd)
  resGDS<-sapply(resL,function(.x)apply(.x$gene_disease_separation[,3:dim(.x$gene_disease_separation)[2]],c(1,2),as.numeric),simplify = "array")
  m<-apply(resGDS,c(1,2),mean,na.rm =TRUE)
  RANds<-cbind(as.data.frame(resL[[1]]$gene_disease_separation[,1:2]),as.data.frame(m))
  ##comment out for moment
  CN <- colnames(ds)[3:length(ds[1,])]

  ##--- output results file comparing observed disease pairs against randomised distribution.
disease_location_sig           <- matrix(0 ,ncol=7, nrow=length(CN))
colnames(disease_location_sig) <- c("HDO.ID","N","mean_ds","SD_ds","Ran_mean_ds","Ran_SD_ds","Utest.pvalue")
disease_location_sig[,1]       <- CN
disease_location_sig[,2]<-loc[match(disease_location_sig[,1],loc[,1]),2]

  ## significance of ds for each disease
  for( i in 1:length(CN) ){

    ## gda matching indices
    indx <- ds[,(2+i)]!="."

    ## gene ids
    ids <- ds[indx,1]

    ## observed ds values
    DS       <- as.numeric(as.vector(ds[indx,(2+i)]))
    disease_location_sig[i,3] <- as.numeric(mean(DS))
    disease_location_sig[i,4] <- as.numeric(sd(DS))

    indy <- match(ids,RANds[,1])

    ## random ds values
    RDS <- as.numeric(as.vector(RANds[indy,(2+i)]))
    disease_location_sig[i,5] <- as.numeric(mean(RDS,na.rm=TRUE))
    disease_location_sig[i,6] <- as.numeric(sd(RDS,na.rm=TRUE))
    disease_location_sig[i,7] <- 1.0

    ## compute wilcox test between observable ds and random ds, and store p.values,
    ## see (Menche et al., 2015).
    if( !is.infinite(DS) && !is.nan(DS) && !is.na(DS) &&  !is.infinite(RDS) && !is.nan(RDS) && !is.na(RDS) ){
      if( length(DS) != 0 && length(RDS) != 0 ){
        wt       <- wilcox.test(DS,RDS)
        disease_location_sig[i,7] <- as.numeric(wt$p.value)
      }
    }
  }

}
