##
# Calculate each diease-pair overlap/seperation on a selected  
# synaptic PPI network models, based on analysis described in: 
# Menche, J. et al. Uncovering disease-disease relationships through the incomplete interactome.
# Science, 347, (6224):1257601 (2015).
##

source('../setUp.R')

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

#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

dadir <- sprintf("%s/%s",subDIR[S],gdaDIR[gdas])
if( !file_test("-d",dadir) ){
    dir.create(dadir)
}
#---



#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

#--- Find all Gene Disease Associations
GDA <- V(gg)$TopOntoOVG

#--- The number of GDA's in graph
NN  <- length(which(GDA!=""))

#---Check
#---Remove Diseases with zero GDA's
remove <- c()
cat("Following Diseases have zero GDA's.\n")
for( d in 1:length(dtype) ){
    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    if( length(IDS) == 0 ){
        cat(dtype[d], " => ", length(IDS),"\n")
	remove <- c(remove,d)
     }
}

if( length(remove) == 0 ){ cat(" => None.\n"); }

if( length(remove) > 0 ){
   disn  <- disn [-remove]	
   dtype <- dtype[-remove]
   disl  <- disl[-remove]
}    
#---     

#--- output result file
res           <- matrix(0 ,ncol=4, nrow=length(dtype))
colnames(res) <- c("Disease","N","mean_ds","SD_ds")
res[,1]       <- dtype


#--- store minimum shorest paths for each gda, and each disease
oo <- matrix(".",nrow=NN,ncol=(length(dtype)+2))
colnames(oo) <- c("Gene.ID","Gene.Name",dtype)
oo[,1]       <- V(gg)$name[GDA !=""]
oo[,2]       <- V(gg)$GeneName[GDA !=""]

##--- loop over each disease
for( d in 1:length(dtype) ){

    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
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

#outfile, the disase localisation info, mean shortest ds for each disease
outfile <- file(sprintf("%s/%s_disease_localisation.csv",dadir,subDIR[S]),"w")
write.table(res, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

#outfile, the gda shortest distance info
outfile <- file(sprintf("%s/%s_gene_disease_separation.csv",dadir,subDIR[S]),"w")
write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

#Disease-disease overlap
DAB <- matrix(".",ncol=length(dtype),nrow=length(dtype))
colnames(DAB) <- dtype
rownames(DAB) <- dtype

#--- NOTE ---#
# DAB is bound by -dmax <= DAB <= dmax
# where dmax denotes the diameter of the network
# dmax <- diameter(gg,directed=F)
#------------#

##--- calculate disease-disease overlap
for( i in 1:length(dtype) ){
    for( j in i:length(dtype) ){

        DAB[i,j] <- 0
        
        if( i != j ){
            DAB[i,j] <- diseaseOverlap(gg,GDA,rownames(DAB)[i],colnames(DAB)[j],oo)
        }
        
    }
}

##--- output file, the disase pair overlap info
outfile <- file(sprintf("%s/%s_disease_separation.csv",dadir,subDIR[S]),"w")
write.table(DAB, file=outfile, append=T, row.names=T, col.names=T, sep="\t",quote=F);
close(outfile);



