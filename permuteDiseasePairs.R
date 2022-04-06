library(igraph);

## Overlap of Disease A and B in the interactome
diseaseOverlap <- function(GG, disA, disB, OO){

## disease A genes 
indA  <- which(colnames(OO)==disA)    
IDS1  <- as.character(OO[OO[,indA[1]]!=".",1])
NIDS1 <- length(IDS1)

## disease B genes 
indB  <- which(colnames(OO)==disB)    
IDS2  <- as.character(OO[OO[,indB[1]]!=".",1])
NIDS2 <- length(IDS2)

## disease A given B
paths  <- igraph::shortest.paths(GG,IDS1,IDS2,weights=NA)
dsA    <- as.numeric(as.vector(apply(paths,1,min))) 

## disease B given A
paths <- igraph::shortest.paths(GG,IDS2,IDS1,weights=NA) 
dsB   <- as.numeric(as.vector(apply(paths,1,min)))	     

## network-based separation between disease A and B 
dAB <- (sum(dsA)+sum(dsB))/(NIDS1+NIDS2)

## network-based localisation of disease A
dA   <- mean(as.numeric(as.vector(OO[OO[,indA[1]]!=".",indA[1]])))

## network-based localisation of disease B
dB   <- mean(as.numeric(as.vector(OO[OO[,indB[1]]!=".",indB[1]])))
    
## overlap between disease A and B
sAB = as.numeric(dAB) - (as.numeric(dA)+as.numeric(dB))/2

return(sAB)
        
}

## permute biological condition data, i.e
## the sample (given by the cols) and genes
## (given by the rows).
## This needs to be run on 'eddie' ECDF, since ~10,000 iterations needed.
permute <- function(GNS, N){

	temp <- sample(GNS,N,replace=F)

	return(temp)
    
}

##---Script required inputs
args  <- commandArgs(TRUE);
SEED  <- as.numeric(args[1]) #random seed no. 
ITS   <- as.numeric(args[2]) #no: of iterations

## set the random number seed
set.seed(as.numeric(SEED))

cat("Running DiseaseLoc.R: \n")
cat("SEED: ", SEED, "\n")
cat("Permutations: ", ITS, "\n")

##---YOUR DATASET TO READ IN
files   <- list.files()
files   <- files[grepl(".gml" ,files,fixed=T)]

studies <- unlist(strsplit(files,".gml"))

##--PPI networks
gg <- read.graph(files[1],format="gml")

## Get all gene Entrez IDS
GNS <- V(gg)$name

##---Load HDO ID DISEASES of INTEREST
##---Use HDO Disease short names
##   Note this information 
source('annotationTYPES.R')

## load gda's for the ppi network
GDA <- V(gg)$TopOntoOVG

##---Remove Diseases with zero GDA's
remove <- c()
for( d in 1:length(dtype) ){
    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    if( length(IDS) == 0 ){
        cat(dtype[d], " => ", length(IDS),"\n")
	remove <- c(remove,d)
     }
}
if( length(remove) > 0 ){
   disn  <- disn [-remove]	
   dtype <- dtype[-remove]
}    
#---  

## for each gda, find the minimum shortest path to next gda (of the same disease)
oo <- matrix(0,nrow=length(GNS),ncol=(length(dtype)+2))
colnames(oo) <- c("Gene.ID","Gene.Name",dtype)
oo[,1] <- V(gg)$name
oo[,2] <- V(gg)$GeneName

N     <- length(dtype)
NELE  <- N*(N+1)/2

## container for mean of Disease-disease overlap
DABm <- matrix(0,nrow=NELE, ncol=3)

## container for sd of Disease-disease overlap
DABsd <- matrix(0,nrow=NELE, ncol=3)

## container for disease-disease overlaps
DDints <- matrix(0, nrow=NELE, ncol=(2+ITS))

for( k in 0:(NELE-1) ){

  ##--- indexing for symmetric matrix
  i = floor( (2*N+1 - sqrt( (2*N+1)*(2*N+1) - 8*k ))/2 );
  j = k - N*i + i*(i-1)/2;

  i = i + 1;
  j = j + i;

  DABm[(k+1),1] <- dtype[i]
  DABm[(k+1),2] <- dtype[j]

  DABsd[(k+1),1] <- dtype[i]
  DABsd[(k+1),2] <- dtype[j]

  DDints[(k+1),1] <- dtype[i]
  DDints[(k+1),2] <- dtype[j]

 }


##start timing
ptm <- proc.time()

cat("Running...\n")
for( p in 1:ITS ){

      cat("Permutation No: ", p, "\n")

    ##--- store minimum shorest paths for each gda, and each disease, for
    ##    permutation 'p'
     temp <- matrix(".",nrow=length(GNS),ncol=(length(dtype)+2))
     colnames(temp) <- c("Gene.ID","Gene.Name",dtype)
     temp[,1] <- V(gg)$name
     temp[,2] <- V(gg)$GeneName

     ##--- loop over each disease
     for( d in 1:length(dtype) ){

         ## get gene-disease associations (GDAs)
    	 IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    	 NN  <- length(IDS)

         ## permute the N GDA's relative to all gene ids
         IDS <- permute(GNS, NN) #case

         ## for each gda, find the minimum shortest path to next gda (of the same 
         XX               = igraph::shortest.paths(gg,IDS,IDS,weights=NA)
         diag(XX)         = NA
         ds               = apply(XX,1,min,na.rm=T)
         indX             = match(names(ds),temp[,1])
         temp[indX,(2+d)] = as.vector(ds)
         oo[indX,(2+d)]   = as.numeric(oo[indX,(2+d)]) + as.numeric(ds)

         rm(IDS,XX)

  }
  

  ##--- calculate disease-disease overlap
  for( k in 1:NELE ){

    overlap = 0;
    if(DDints[k,1] != DDints[k,2]){
       overlap <- as.numeric(diseaseOverlap(gg,DDints[k,1],DDints[k,2],temp))
    }

    DDints[k,(2+p)] <- overlap

  }#done 

  }#permutations

pet <- proc.time() - ptm
cat("Finished! ", sprintf("time = %.3f \n", pet[[1]]), "\n")

## calculate the mean and sd of each disease pair 
 for( k in 1:NELE ){

     if(DDints[k,1] != DDints[k,2]){  
        DABm[k,3]  <- as.character(mean(as.numeric(DDints[k,3:(2+ITS)])))
        DABsd[k,3] <- as.character(sd(as.numeric(DDints[k,3:(2+ITS)])))
        } else {
	  DABm[k,3]  <- 0
          DABsd[k,3] <- 0
	}

}#done

outfile <- file("gene_disease_separation.csv","w")
write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

##--- output file, each disease pair overlap for each permutation
outfile <- file("sAB_random_separation.csv","w")
write.table(DDints, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

##--- output file, mean of each disease pair overlap
outfile <- file("mean_disease_separation.csv","w")
write.table(DABm, file=outfile, append=T, row.names=T, col.names=T, sep="\t",quote=F);
close(outfile);

##--- output file, sd of each disease pair overlap
outfile <- file("sd_disease_separation.csv","w")
write.table(DABsd, file=outfile, append=T, row.names=T, col.names=T, sep="\t",quote=F);
close(outfile);



