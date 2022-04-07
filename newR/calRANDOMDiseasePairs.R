##
# Calculate each diease-pair overlap/seperation on a selected  
# synaptic PPI network models, based on analysis described in: 
# Menche, J. et al. Uncovering disease-disease relationships through the incomplete interactome.
# Science, 347, (6224):1257601 (2015).
##

source('../setUp.R')
library(WGCNA) ## for it's qvalue calculation

qscore <- function(zz,FDR){

    LL <- FDR[FDR[,1] < as.numeric(zz),2]

    if( length(LL) != 0 ){ return(LL[end(LL)[1]]); }

    return(1)
}

##---Check or create output dir called RESULTS
if( !file_test("-d","RESULTS") ){
    dir.create("RESULTS")
}

resdir1 <- sprintf("RESULTS/%s",subDIR[S])

if( !file_test("-d",resdir1) ){
    dir.create(resdir1)
}

gdaDIR <- "ovg"
resdir <- sprintf("%s/%s",resdir1,gdaDIR)

if( !file_test("-d",resdir) ){
    dir.create(resdir)
}
#---

FILES <- vector(length=4)
FILES[1] <- "gene_disease_separation"  
FILES[2] <- "disease_separation"
FILES[3] <- "random_separation"
FILES[4] <- "disease_localisation"

## For Disease Sizes on the interactome
loc   <- read.table(sprintf("%s/%s/%s_%s.csv",subDIR[S],gdaDIR,subDIR[S],FILES[4]),sep="\t",header=T)

ds    <- read.table(sprintf("%s/%s/%s_%s.csv",subDIR[S],gdaDIR,subDIR[S],FILES[1]),sep="\t",header=T)

## load randomised gene disease separation data
RANds <- read.table(sprintf("%s/Random_Disease_Pairs/%s/random_%s.csv",rndDIR[1],subDIR[S],FILES[1]),sep="\t",header=T)

##comment out for moment
CN <- colnames(ds)[3:length(ds[1,])]

##--- output results file comparing observed disease pairs against randomised distribution.
res           <- matrix(0 ,ncol=9, nrow=length(CN))
colnames(res) <- c("HDO.ID","Disease.long","Disease","N","mean_ds","SD_ds","Ran_mean_ds","Ran_SD_ds","Utest.pvalue")

res[,3]       <- CN
res[,1]       <- disn[match(res[,3],dtype)]
res[,2]       <- disl[match(res[,3],dtype)]
res[,4]       <- loc[match(res[,3],loc[,1]),2]


## significance of ds for each disease
for( i in 1:length(CN) ){
     
     ## gda matching indices
     indx <- ds[,(2+i)]!="."
     
     ## gene ids
     ids <- ds[indx,1]

     ## observed ds values
     DS       <- as.numeric(as.vector(ds[indx,(2+i)]))
     res[i,5] <- as.numeric(mean(DS))
     res[i,6] <- as.numeric(sd(DS))

     indy <- match(ids,RANds[,1])
     
     ## random ds values
     RDS <- as.numeric(as.vector(RANds[indy,(2+i)]))
     res[i,7] <- as.numeric(mean(RDS))
     res[i,8] <- as.numeric(sd(RDS))
     res[i,9] <- 1.0
    
    ## compute wilcox test between observable ds and random ds, and store p.values, 
    ## see (Menche et al., 2015).
    if( !is.infinite(DS) && !is.nan(DS) && !is.na(DS) &&  !is.infinite(RDS) && !is.nan(RDS) && !is.na(RDS) ){
        if( length(DS) != 0 && length(RDS) != 0 ){ 
            wt       <- wilcox.test(DS,RDS)
            res[i,9] <- as.numeric(wt$p.value)
        }
    }
}

outfile <- file(sprintf("%s/Disease_location_sig.csv",resdir),"w");
write.table( res, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);


##--- load observed mean sAB, i.e. overlap between disease A and B
sAB <- read.table(sprintf("%s/%s/%s_%s.csv",subDIR[S],gdaDIR,subDIR[S],FILES[2]),sep="\t",header=T)

##--- load random sAB values for each iteration
RAW_sAB  <- read.table(sprintf("%s/Random_Disease_Pairs/%s/RAW_random_sAB_%s.csv",rndDIR[1],subDIR[S],FILES[3]),sep="\t",header=T)

##---no. of random permutations
perms <- length(RAW_sAB[1,])

##--- save random mean and sd sAB values in containers 
RAN_sAB_mean <- matrix(0,ncol=3,nrow=length(RAW_sAB[,1]))
RAN_sAB_sd   <- matrix(0,ncol=3,nrow=length(RAW_sAB[,1]))

RAN_sAB_mean[,1] <- as.character(RAW_sAB[,1])
RAN_sAB_mean[,2] <- as.character(RAW_sAB[,2])

RAN_sAB_sd[,1] <- as.character(RAW_sAB[,1])
RAN_sAB_sd[,2] <- as.character(RAW_sAB[,2])
##---

for( i in 1:length(RAW_sAB[,1]) ){
    RAN_sAB_mean[i,3] <- as.character(mean(as.vector(unlist(RAW_sAB[i,3:perms]))))
    RAN_sAB_sd[i,3]   <- as.character(sd(as.vector(unlist(RAW_sAB[i,3:perms]))))
}

Nn    <- length(CN)
NELE  <- Nn*(Nn+1)/2

##---no: of levels for Bonferroni correction
Nlevels = NELE;

##--- Output file for disease-disease separation/overlap
CN <-  c("HDO.ID","Disease.long","Disease","N","HDO.ID","Disease.long","Disease","N","sAB","Separated","Overlap","zScore","pvalue","Separation/Overlap.than.chance","Bonferroni","p.adjusted","q-value")
zs <- matrix(".", nrow=NELE, ncol=length(CN))
colnames(zs) <- CN

#
for( k in 0:(NELE-1) ){

  ##--- linear indexing for symmetric matrix
  i = floor( (2*Nn+1 - sqrt( (2*Nn+1)*(2*Nn+1) - 8*k ))/2 );
  j = k - Nn*i + i*(i-1)/2;

  i = i + 1;
  j = j + i;

  zScore = 0

    if( !is.nan(as.numeric(RAN_sAB_sd[(k+1),3])) ){
    
  ## compute z-score, i.e. separation of mean sAB, against a randomised model (of the mean of sAB),   
  ## see (Menche et al., 2015).
  if( as.numeric(RAN_sAB_sd[(k+1),3]) != 0){
      zScore = (as.numeric(as.vector(sAB[i,j])) - as.numeric(as.vector(RAN_sAB_mean[(k+1),3])))/(as.numeric(as.vector(RAN_sAB_sd[(k+1),3])))
  }

  ## compute p.value from the normal distribution
  ## See also http://www.cyclismo.org/tutorial/R/pValues.html     
  pval <- pnorm(-abs(zScore))
  pval <- 2 * pval

  zs[(k+1),1] <- as.character(disn[which(dtype==RAN_sAB_mean[(k+1),1])])
  zs[(k+1),2] <- as.character(disl[which(dtype==RAN_sAB_mean[(k+1),1])])    
  zs[(k+1),3] <- as.character(RAN_sAB_mean[(k+1),1])
  zs[(k+1),4] <- as.character(loc[which(loc[,1]==RAN_sAB_mean[(k+1),1]),2])

  zs[(k+1),5] <- as.character(disn[which(dtype==RAN_sAB_mean[(k+1),2])])
  zs[(k+1),6] <- as.character(disl[which(dtype==RAN_sAB_mean[(k+1),2])])    
  zs[(k+1),7] <- as.character(RAN_sAB_mean[(k+1),2])
  zs[(k+1),8] <- as.character(loc[which(loc[,1]==RAN_sAB_mean[(k+1),2]),2])
    
  ## sAB, the disease-disease separation/overlap measure, on the interactome   
  zs[(k+1),9] <- as.character(sAB[i,j])

  ## sAB > 0, implies separation
  zs[(k+1),10] <- ifelse((as.numeric(zs[(k+1),9]) > 0), "YES", ".")

  ## sAB < 0, implies overlap
  zs[(k+1),11] <- ifelse((as.numeric(zs[(k+1),9]) < 0), "YES", ".")

  ## save z-score and p.value
  zs[(k+1),12] <- as.character(zScore)
  zs[(k+1),13] <- as.character(pval)

  ## z-scores < 0 (>0), implies separation/overlap smaller (larger) than by chance  
  zs[(k+1),14] <- ifelse((as.numeric(zs[(k+1),12]) < 0), "Smaller", "larger")

  ## Bonferroni correction for p.value ('stars' can be found in 'setUp.R')
  temp <- "."
  for( x in 1:length(alpha) ){
      if(as.numeric(zs[(k+1),13]) < as.numeric(alpha[x]/Nlevels)){ temp <- stars[x] }
  }

  ## save the Bonerroni correction, repersented by stars, here.      
  zs[(k+1),15] <- temp

    ## default fill of output container 
    } else {
        zs[(k+1),1] <- as.character(disn[which(dtype==RAN_sAB_mean[(k+1),1])])
        zs[(k+1),2] <- as.character(disl[which(dtype==RAN_sAB_mean[(k+1),1])])    
        zs[(k+1),3] <- as.character(RAN_sAB_mean[(k+1),1])
        zs[(k+1),4] <- as.character(loc[which(loc[,1]==RAN_sAB_mean[(k+1),1]),2])
        
        zs[(k+1),5] <- as.character(disn[which(dtype==RAN_sAB_mean[(k+1),2])])
        zs[(k+1),6] <- as.character(disl[which(dtype==RAN_sAB_mean[(k+1),2])])    
        zs[(k+1),7] <- as.character(RAN_sAB_mean[(k+1),2])
        zs[(k+1),8] <- as.character(loc[which(loc[,1]==RAN_sAB_mean[(k+1),2]),2]) 
    }

}
  
##--- pre-calculate each random z-score, using standard deviation (sd) obtained from each eddie job...
##    you could use here the global sd obtained from all random permutations.

##    container to save randomised z-scores
tests <- matrix(0, nrow=length(RAW_sAB[,1]),ncol=(length(RAW_sAB[1,])-2))
for( k in 0:(NELE-1) ){

    ##--- the linerar indexing for symmetric matrix
    i = floor( (2*Nn+1 - sqrt( (2*Nn+1)*(2*Nn+1) - 8*k ))/2 );
    j = k - Nn*i + i*(i-1)/2;
    
    i = i + 1;
    j = j + i;

    ##
    if( zs[(k+1),1] != zs[(k+1),4] ){
  
        test <- 0
          
        Mn   <- as.numeric(RAN_sAB_mean[(k+1),3])
        Sd   <- as.numeric(RAN_sAB_sd[(k+1),3])
        
        RsAB <- as.numeric(RAW_sAB[(k+1),-c(1,2)])
            
        ## store the random z-score 
        tests[(k+1),] <- (as.numeric(RsAB - Mn))/Sd
        
    }

}

##--- save random zscores for plotting    
outfile <- file(sprintf("%s/random_zscores.csv",resdir),"w");
write.table( tests, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);



##--- FDR calculation ---##
# Note: In my study i've only 66 disease-disease pairs to consider when correcting for multiple hypothesis
#       corrections, and i focus on Bonerroni corrections @ 0.001, implying i expect to get <<1 false positive.
#       If you've >> 100 pairs to consider, you'll probably need to use the FDR calculation below. 
##-----------------------## 

## switch 'calFDR', should we calculate FDR?
calFDR=0
##calFDR=1

## save p.adjusted value in output container
zs[,16] <- p.adjust(as.numeric(zs[,13]),method="BY")

## if calFDR is FALSE, we'll use WGCNA's qvalue calculation for FDR
zs[,17] <- qvalue(as.numeric(zs[,13]))$qvalue

if(calFDR){
    
    ##--- Calculate local and Global FDR 
    zDeltas    <- seq(0,10,0.01)    
    FDrate     <- matrix(0,ncol=2,nrow=length(zDeltas))
    FDrate[,1] <- zDeltas

    ##--- loop through each zDelta and test the fraction of random z-scores > than it,
    ##    this will give use the local FDR for the range (i.e. zDeltas) of z-scores, which
    ##    we can then compare with the observed z-scores
    for( f in 1:length(FDrate[,1]) ){
        
        zthres <- rep(0,length(zs[,1]))

        ## norm for global FDR
        norm <- 0

        for( k in 0:(NELE-1) ){

            ##--- linear indexing for symmetric matrix
            i = floor( (2*Nn+1 - sqrt( (2*Nn+1)*(2*Nn+1) - 8*k ))/2 );
            j = k - Nn*i + i*(i-1)/2;
            
            i = i + 1;
            j = j + i;

            if( zs[(k+1),1] != zs[(k+1),5] ){

     
                norm = norm + 1

                ## number of random z-scores > than the z-score under investigation
                tally = 0
                tally = length(which((abs(tests[(k+1),]) > FDrate[f,1]) == TRUE))

                ##---local FDR at zthreshold
                zthres[(k+1)] <- as.numeric(tally)/perms
      
            }

        }

        ##--- Global FDR at zthreshold = zDeltas[f] 
        FDrate[f,2] <- (1/norm)*(sum(as.numeric(zthres)))

    }

    ##--- At this point you can use FDrate to calculate the qvalue, the minimal global FDR at which a disease-disease pair becomes
    ##    reaches significances.    
    Zz      <- abs(as.numeric(zs[,12]))
    qs      <- lapply(Zz,qscore,FDR=FDrate)
    qs      <- unlist(qs)
    zs[,16] <- qs
    
}#if

## save output file
outfile <- file(sprintf("%s/Disease_overlap_sig.csv",resdir),"w");
write.table( zs, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);
