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

