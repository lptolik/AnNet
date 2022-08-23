
source('../setUp.R');

##--- Load clusterCons library
library(clusterCons)

##--- Get cluster robustness values (usig R's clusterCons package)
##--- run clusterCons's patch called 'memrob.R' (this script is needed)
source('memrob.R');

##--- Directories needed
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]

##--- Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

##--- declare clustering algorithms used
##--- Clustering algorithms used
alg <- "Spectral"

##--- load corresponding graph which was used to build the consensus matrices from
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

##--- build reference matrix from the graph
refin     <- matrix("",ncol=3,nrow=length(V(gg)$name))
refin[,1] <- igraph::get.vertex.attribute(gg,"name",V(gg))
refin[,2] <- igraph::get.vertex.attribute(gg,"name",V(gg))
refin[,3] <- igraph::get.vertex.attribute(gg,alg,V(gg))

##--- path to consensus matrix
st1 = sprintf("%s/%s/%s/consensusmatrix.txt.gz",rndDIR[1],subDIR[S],alg[a]);

##--- Read in consensus matrix
filein = read.table(gzfile(st1), header=FALSE, sep=",");
dimnames(filein)[2] <- dimnames(filein)[1]

##--- format reference matrix
##--- the reference matrix with the correct row.names(as a data.frame)
refmat       <- as.matrix(refin);
ref          <- refmat[,2:3];
rm           <- data.frame(ref);
rownames(rm) <- rm$X1;
rm$X1        <- NULL;
names(rm)    <- 'cm';

##--- format consensus matrix
##--- the consensus matrix you may have made (as a numeric matrix)
conmat       <- as.matrix(filein);
cm           <- data.frame(conmat);
names(cm)    <- rownames(rm);
rownames(cm) <- rownames(rm);
cm           <- as.matrix(cm);

##--- max number of clusters
kk     = max(as.numeric(as.vector(rm$cm)));

##--- make the consensus matrix object for clusterCons so you can use its functions
out <- new('consmatrix', cm=cm, rm=rm, k=kk, a=alg);

##--- get cluster robustness values from clusterCons
cr <- clusterCons::clrob(out);

##--- the scaled cluster robustness values
crScales <- cr$rob
crScales <- (crScales-min(crScales))/(max(crScales)-min(crScales))

##--- create output file
oo <- data.frame(a=as.character(),b=as.numeric(),c=as.numeric(),d=as.numeric(), e=as.numeric())
oo <- rbind(oo,data.frame(a=as.character(rep(alg,length(rownames(cr)))),
                          b=as.numeric(rownames(cr)),
                          c=as.numeric(table(rm[,1])),
                          d=as.numeric(cr$rob),
                          e=as.numeric(crScales)))

##--- check subDir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

##--- save output file
colnames(oo) <- c("alg","C","Cn","Crob","CrobScaled")
outfile <- file(sprintf("%s/%s_ClusterRobustness.csv",subDIR[S],subDIR[S]),"w")
write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);




