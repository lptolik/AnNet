#' @export
#' @import WGCNA
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

stars    <- c("*","**","***")

#' Auxiliary function to replase NAs with zeros.
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
zeroNA<-function(x){
  x[is.na(x)]<-0
  return(x)
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
    map=cbind(map,ifelse(grepl(dtype[i],GDA,fixed=TRUE),1,0))
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
  if(permute!='none'){
    rgda <- matrix(NA,nrow=vcount(gg),ncol=(length(diseases)))
    colnames(rgda) <- c(diseases)
    if(permute=='binned'){
      map <- degree.binned.GDAs(gg,gda,diseases)
    }

    for( d in 1:length(diseases) ){

      IDS <- V(gg)$name[grepl(diseases[d],gda,fixed=T)]
      N   <- length(IDS)

      if(permute=='random'){
        ## permute the N GDA's relative to all gene ids
        IDS <- AnNet::permute( V(gg)$name, N) #case
      }else if(permute=='binned'){
        IDS <- sample.deg.binned.GDA(map,diseases[d])
      }
      rgda[match(IDS,V(gg)$name),d]<-1
    }
    gda<-apply(rgda,1,function(.x)paste(diseases[!is.na(.x)],collapse = COLLAPSE))
  }
  res           <- matrix(0 ,ncol=4, nrow=length(diseases))
  colnames(res) <- c("Disease","N","mean_ds","SD_ds")
  res[,1]       <- unescapeAnnotation(diseases)


  #--- store minimum shorest paths for each gda, and each disease
  oo <- matrix(".",nrow=length(gda),ncol=(length(diseases)+2))
  colnames(oo) <- c("Gene.ID","Gene.Name",diseases)
  oo[,1]       <- V(gg)$name#[gda !=""]
  oo[,2]       <- V(gg)$GeneName#[gda !=""]

  ##--- loop over each disease
  for( d in 1:length(diseases) ){

    IDS <- V(gg)$name[grepl(diseases[d],gda,fixed=T)]
    N   <- length(IDS)

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
  colnames(DAB) <- unescapeAnnotation(diseases)
  rownames(DAB) <- unescapeAnnotation(diseases)
  colnames(oo) <- c("Gene.ID","Gene.Name",unescapeAnnotation(diseases))

  if(permute=='none'){
    oo<-oo[gda !="",]
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
runPermDisease<-function(gg,name,diseases=NULL,Nperm=100,alpha=c(0.05,0.01,0.001)){
  mean0<-function(x){
    return(mean(zeroNA(x)))
  }
  sd0<-function(x){
    return(sd(zeroNA(x)))
  }
  resD<-calcDiseasePairs(gg=gg,name=name,diseases=diseases,permute = 'none')
  ds<-resD$gene_disease_separation
  loc<-resD$disease_localisation
  resL<-lapply(1:Nperm,function(.x)calcDiseasePairs(gg=gg,name=name,diseases=diseases,permute='random'))
  resGDS<-sapply(resL,function(.x)apply(.x$gene_disease_separation[,3:dim(.x$gene_disease_separation)[2]],c(1,2),as.numeric),simplify = "array")
  m<-apply(resGDS,c(1,2),mean0)
  RANds<-cbind(as.data.frame(resL[[1]]$gene_disease_separation[,1:2]),as.data.frame(m))
  ##comment out for moment
  disn <- colnames(ds)[3:length(ds[1,])]

  ##--- output results file comparing observed disease pairs against randomised distribution.
  disease_location_sig           <- matrix(0 ,ncol=7, nrow=length(disn))
  colnames(disease_location_sig) <- c("HDO.ID","N","mean_ds","SD_ds","Ran_mean_ds","Ran_SD_ds","Utest.pvalue")
  disease_location_sig[,1]       <- disn
  disease_location_sig[,2]<-loc[match(disease_location_sig[,1],loc[,1]),2]

  ## significance of ds for each disease
  for( i in 1:length(disn) ){

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
    disease_location_sig[i,5] <- as.numeric(mean0(RDS))
    disease_location_sig[i,6] <- as.numeric(sd0(RDS))
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
  sAB<-resD$disease_separation
  RAW_sAB<-sapply(resL,function(.x).x$disease_separation,simplify = "array")
  RAN_sAB_mean<-apply(RAW_sAB,c(1,2),mean0)
  RAN_sAB_sd<-apply(RAW_sAB,c(1,2),sd0)
  perms <- dim(RAW_sAB)[3]
  Nn    <- length(disn)
  NELE  <- Nn*(Nn+1)/2

  ##---no: of levels for Bonferroni correction
  Nlevels = NELE;

  ##--- Output file for disease-disease separation/overlap
  #CN <-  c("HDO.ID","Disease.long","Disease","N","HDO.ID","Disease.long","Disease","N","sAB","Separated","Overlap","zScore","pvalue","Separation/Overlap.than.chance","Bonferroni","p.adjusted","q-value")
  CN <-  c("HDO.ID","N","HDO.ID","N","sAB","Separated","Overlap","zScore","pvalue","Separation/Overlap.than.chance","Bonferroni","p.adjusted","q-value")
  zs <- matrix(".", nrow=NELE, ncol=length(CN))
  colnames(zs) <- CN
  tests <- matrix(0, nrow=NELE,ncol=perms)
  for( k in 0:(NELE-1) ){

    ##--- linear indexing for symmetric matrix
    i = floor( (2*Nn+1 - sqrt( (2*Nn+1)*(2*Nn+1) - 8*k ))/2 );
    j = k - Nn*i + i*(i-1)/2;

    i = i + 1;
    j = j + i;

    zScore = 0

    if( !is.nan(as.numeric(RAN_sAB_sd[i,j])) ){

      ## compute z-score, i.e. separation of mean sAB, against a randomised model (of the mean of sAB),
      ## see (Menche et al., 2015).
      if( as.numeric(RAN_sAB_sd[i,j]) != 0){
        zScore = (as.numeric(as.vector(sAB[i,j])) - as.numeric(as.vector(RAN_sAB_mean[i,j])))/(as.numeric(as.vector(RAN_sAB_sd[i,j])))
      }

      ## compute p.value from the normal distribution
      ## See also http://www.cyclismo.org/tutorial/R/pValues.html
      pval <- pnorm(-abs(zScore))
      pval <- 2 * pval

      zs[(k+1),1] <- disn[i]
      zs[(k+1),2] <- as.character(loc[which(loc[,1]==disn[i]),2])

      zs[(k+1),3] <- disn[j]
      zs[(k+1),4] <- as.character(loc[which(loc[,1]==disn[j]),2])

      ## sAB, the disease-disease separation/overlap measure, on the interactome
      zs[(k+1),5] <- as.character(sAB[i,j])

      ## sAB > 0, implies separation
      zs[(k+1),6] <- ifelse((as.numeric(zs[(k+1),5]) > 0), "YES", ".")

      ## sAB < 0, implies overlap
      zs[(k+1),7] <- ifelse((as.numeric(zs[(k+1),5]) < 0), "YES", ".")

      ## save z-score and p.value
      zs[(k+1),8] <- as.character(zScore)
      zs[(k+1),9] <- as.character(pval)

      ## z-scores < 0 (>0), implies separation/overlap smaller (larger) than by chance
      zs[(k+1),10] <- ifelse((as.numeric(zs[(k+1),8]) < 0), "Smaller", "larger")

      ## Bonferroni correction for p.value ('stars' can be found in 'setUp.R')
      temp <- "."
      for( x in 1:length(alpha) ){
        if(as.numeric(zs[(k+1),9]) < as.numeric(alpha[x]/Nlevels)){ temp <- stars[x] }
      }

      ## save the Bonerroni correction, repersented by stars, here.
      zs[(k+1),11] <- temp
      ## default fill of output container
    } else {
      zs[(k+1),1] <- disn[i]
      zs[(k+1),2] <- as.character(loc[which(loc[,1]==disn[i]),2])

      zs[(k+1),3] <- disn[j]
      zs[(k+1),4] <- as.character(loc[which(loc[,1]==disn[j]),2])
    }
    # if( zs[(k+1),1] != zs[(k+1),3] ){
    #
    #   test <- 0
    #
    #   Mn   <- as.numeric(RAN_sAB_mean[i,j])
    #   Sd   <- as.numeric(RAN_sAB_sd[i,j])
    #
    #   RsAB <- as.numeric(RAW_sAB[i,j,])
    #
    #   ## store the random z-score
    #   tests[(k+1),] <- (as.numeric(RsAB - Mn))/Sd
    #
    # }

  }
  ## save p.adjusted value in output container
  zs[,12] <- p.adjust(as.numeric(zs[,9]),method="BY")

  ## if calFDR is FALSE, we'll use WGCNA's qvalue calculation for FDR
  zs[,13] <- qvalue(as.numeric(zs[,9]))$qvalue

  return(list(Disease_overlap_sig=zs,Disease_location_sig=disease_location_sig))
}
