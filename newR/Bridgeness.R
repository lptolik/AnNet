source('../setUp.R');


#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("POWERlawFIT",DIRS)]
OUT[3] <- DIRS[grepl("SVI",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}


#---Set Options
runBridge  <- vector(length=2)
runBridge[1] <- 1  #Calculate Bridgeness
runBridge[2] <- 0  #Plot Bridgeness

#---declare clustering algorithms in graph, and with a corresponding consensus matrix
alg <- c("Spectral")

#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

if(runBridge[1]){

    #---number of vertices/genes
    N    <- length(V(gg))

    #---number of edges/PPIs
    M    <- length(E(gg))
    
    #---container column names
    CN   <- c('ENTREZ.ID','GENE.NAME',sprintf("BRIDGENESS.%s",alg))

    #---container to store Bridgeness for algorithm 'alg'
    meas <- matrix(0, nrow=N, ncol=length(CN))
    colnames(meas) <- CN

    meas[,1] <- as.character(V(gg)$name)
    meas[,2] <- as.character(V(gg)$GeneName)
      
    ##START filling meas after PageRank column
    FROM <- 2

    ## run over each algorithm, only one in this example
    for( a in 1:length(alg) ){

        cat("calculating Bridgeness for: ", alg[a], "\n")
    
        ##--- build reference matrix from the graph
        refin     <- matrix("",ncol=3,nrow=length(V(gg)$name))
        refin[,1] <- igraph::get.vertex.attribute(gg,"name",V(gg))
        refin[,2] <- igraph::get.vertex.attribute(gg,"name",V(gg))
        refin[,3] <- igraph::get.vertex.attribute(gg,alg[a],V(gg))
    
        ##--- path to consensus matrix
        st1 = sprintf("%s/%s/%s/consensusmatrix.txt.gz",rndDIR[1],subDIR[S],alg[a]);
        
        ##--- Read in consensus matrix
        filein = read.table(gzfile(st1), header=FALSE, sep=",");
        dimnames(filein)[2] <- dimnames(filein)[1]
        
        ##format reference matrix
        ##the reference matrix with the correct row.names(as a data.frame)
        refmat = as.matrix(refin);
        ref    = refmat[,2:3];
        rm           <- data.frame(ref);
        rownames(rm) <- rm$X1;
        rm$X1        <- NULL;
        names(rm)    <- 'cm';

        ##format consensus matrix
        ##the consensus matrix you may have made (as a numeric matrix) 
        conmat = as.matrix(filein);
        cm           <- data.frame(conmat);
        names(cm)    <- rownames(rm);
        rownames(cm) <- rownames(rm);
        cm           <- as.matrix(cm);

        rm(filein,refin)

        ##get consensus matrix indices for each edge in edgelist
        indA <- match(igraph::get.edgelist(gg)[,1],rownames(cm))
        indB <- match(igraph::get.edgelist(gg)[,2],rownames(cm))

        dat  <- data.frame(indA,indB)
    
        ##get community assigned to each vertex in edgelist from the algorithm 'alg'
        elA    <- igraph::get.vertex.attribute(gg,alg[a],V(gg))[match(igraph::get.edgelist(gg)[,1],V(gg)$name)]
        elB    <- igraph::get.vertex.attribute(gg,alg[a],V(gg))[match(igraph::get.edgelist(gg)[,2],V(gg)$name)]

        ##for each edge record the community assigned to each vertex and it's consensus matrix value 
        ed      <- matrix(ncol=6,nrow=length(E(gg)))
        ed[,1]  <- igraph::get.edgelist(gg)[,1]
        ed[,2]  <- igraph::get.edgelist(gg)[,2]
        ed[,3]  <- elA
        ed[,4]  <- elB
        ed[,5]  <- apply(dat,1,function(x,mat) mat[x[1],x[2]], mat=cm)
        ed[,6]  <- (elA-elB)

        ##maximum number of communities found by clustering algorithm
        Cmax  <- max(igraph::get.vertex.attribute(gg,alg[a],V(gg)))
    
        rm(cm)
            
        ##loop over each vertex in the graph
        for( i in 1:length(V(gg)) ){

            ##get edges belonging to the i'th veretx
            ind <- which(ed[,1] == V(gg)$name[i] | ed[,2] == V(gg)$name[i])

            ##get community belonging to the i'th vertex       
            c <- igraph::get.vertex.attribute(gg,alg[a],V(gg))[i]

            ##reorder edge communities, so ed[,3] equals current community no: 'c'
            for( k in 1:length(ind) ){
                if( ed[ind[k],6] != 0 && ed[ind[k],4] == c ){
                    ed[ind[k],4] <- ed[ind[k],3]
                    ed[ind[k],3] <- c
                }
            }        
        
            ##number of communities i'th vertex is connected too (via it's edges)
            cc <- unique(ed[ind,4])

            ##use sum of consensus values to calculate the likelihood of i'th
            ##vertex beloning to to k'th community. 
            prob <- vector(length=length(cc))
            for( k in 1:length(cc) ){                
                prob[k] = sum(as.numeric(ed[which(ed[ind,4]==cc[k]),5]))/length(ind)
            }

            ##normalise
            prob <- prob/sum(prob)
        
            ##calculate bridgeness of i'th vertex
            ##Fuzzy communities and the concept of bridgeness in complex networks, T. Nepusz, arXiv, 2007
            b    <- sum( (prob - 1/Cmax) * (prob - 1/Cmax))

            Kzero <- Cmax - length(cc)
            b = b + sum(rep((1/(Cmax*Cmax)),times=Kzero))
        
            ##store values
            ##BRIDGENESS.
            meas[i,(FROM+a)]  <- 1-sqrt( Cmax/(Cmax-1) * b )

        }            

    }
    
    outfile <- file(sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]),"w")
    write.table(meas, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
    close(outfile);

}


    
