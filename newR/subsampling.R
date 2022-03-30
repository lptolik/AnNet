library(igraph);
library(CDMSuite);

#---Find Largest CC
findLCC <- function(GG){

    dec <- igraph::decompose.graph(GG)
    d=1
    CC=length(V(dec[[1]]))
    for( i in 1:length(dec) ){
        if(length(V(dec[[i]])) > CC){
            d=i
            CC=length(V(dec[[i]]))
        }
    }

    GG  <- igraph::decompose.graph(GG)[[d]]
    return(GG)

}

#---Get all edges internal to a community
intraEdges <- function(GG, ALG, CC, INTRA=NULL, INTER=NULL){

    intra = NULL #edges in the community CC
    inter = NULL #edges going out from community CC

    if( !is.null(igraph::get.vertex.attribute(GG,ALG)) ){

        coms <- get.vertex.attribute(GG,ALG)

        if( length(which(coms == CC)) != 0 ){

            ed_cc = E(GG)[inc(coms == CC)]

            all_edges_m <- get.edges(GG, ed_cc) #matrix representation

            inter = (ed_cc[!(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])

            intra = (ed_cc[(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])

        }
    }

    if( INTRA==TRUE && !is.null(intra) ){
        intra_m = get.edges(GG,intra)
        intra   = cbind(V(GG)$name[intra_m[,1]],V(GG)$name[intra_m[,2]])
        return(intra)
    }

    if( INTER==TRUE && !is.null(inter) ){
        inter_m = get.edges(GG,inter)
        inter   = cbind(V(GG)$name[inter_m[,1]],V(GG)$name[inter_m[,2]])
        return(inter)
    }

    return(NULL)

}


recluster <- function( GG, ALGN, CnMAX, CnMIN=1 ){

    if( !is.null(igraph::get.vertex.attribute(GG,ALGN)) ){

    cat("running reclustering...")

        #--- algorithm clustering 1
        ALG1 <- get.vertex.attribute(GG,ALGN,V(GG))
        ALG1 <- cbind(V(GG)$name, ALG1)

        Cn <- table(as.numeric(ALG1[,2]))
        cc <- names(Cn)[Cn > CnMAX]


        RES <- list()
        k=1
        for( i in 1:length(cc) ){

            edCC = intraEdges(GG, ALGN, cc[i], INTRA=TRUE)

            if( !is.null(edCC) ){

                ggLCC    <- graph_from_data_frame(d=edCC, directed=F)

                #fc
                if( ALGN == "fc" ){
                    res      <- igraph::fastgreedy.community(ggLCC)
                    oo       <- cbind(res$names, res$membership)
                }

                #lec
                if( ALGN == "lec" ){
                    lec     <- igraph::leading.eigenvector.community(ggLCC)
                    res     <- igraph::leading.eigenvector.community(ggLCC, start=membership(lec))
                    oo      <- cbind(res$names, res$membership)
                }

                #wt
                if( ALGN == "wt" ){
                    res      <- igraph::walktrap.community(ggLCC)
                    oo       <- cbind(res$names, res$membership)
                }


                #louvain
                if( ALGN == "louvain" ) {
                    res      <- igraph::cluster_louvain(ggLCC)
                    oo       <- cbind(res$names, res$membership)
                }

                #infomap
                if( ALGN == "infomap" ){
                    res      <- igraph::cluster_infomap(ggLCC)
                    oo       <- cbind(res$names, res$membership)
                }

                #sgG1
                if( ALGN == "sgG1" ){
                    res      <- igraph::spinglass.community(findLCC(ggLCC), spins=as.numeric(500),gamma=1)
                    oo       <- cbind(res$names, res$membership)
                }


                #Spectral
                if( ALGN == "Spectral" ){
                    res <- spectral(DF=edCC, CnMIN=CnMIN)
                    oo  <- cbind(res$ID, res$K)
                }

                RES[[k]]      <- oo
                names(RES)[k] <- cc[i]
                k=k+1
            }

        }#for


        if( length(RES) == 0 ){ return(NULL) }


        #--- algorithm clustering 2
        ALG2     <- cbind(ALG1, rep(-1, length(ALG1[,1])))
        indx     <- match(ALG2[,2],cc)
        indx     <- ifelse(is.na(indx),TRUE, FALSE)
        ALG2[,3] <- ifelse(indx, ALG2[,2], ALG2[,3])

        CCmax = max(as.numeric(ALG2[,3]))

        for( i in 1:length(cc) ){

            temp     <- RES[[i]]
            temp[,2] <- as.numeric(temp[,2]) + CCmax

            indx <- match(ALG2[,1],temp[,1])
            indx <- temp[indx,2]

            ALG2[,3] = ifelse(is.na(indx),ALG2[,3],indx)

            CCmax = max(as.numeric(ALG2[,3]))

        }

        #---reorder ALG2[,3]
        N = length(V(GG));

        temp    <- rep(-1, N)
        counter <- min(as.numeric(ALG2[,3]))
        Knew    <- 1;
        Kmax    <- max(as.numeric(ALG2[,3]))

        while( counter <= Kmax ){

            found=FALSE;

            for(v in 1:N ){
                if( as.numeric(ALG2[v,3]) == counter ){
                    temp[v] = Knew;
                    found=TRUE;
                }
            }

            if(found) Knew=Knew+1;

            counter=counter+1;
        }

	cat("...done.\n")

        #---final
        ALG3 <- cbind(ALG2, temp)
        return(ALG3)
    }

    return(NULL)

}

#---Script required inputs
args  <- commandArgs(TRUE);
aa    <- as.numeric(args[1]) #clustering alg.
type  <- as.numeric(args[2]) #edges=>1 or nodes=>2  to mask
mask  <- as.numeric(args[3]) #number of edges/nodes to mask
Cnmin <- as.numeric(args[4]) #Cn min for Spectral algorithm
Cnmax <- as.numeric(args[5]) #Cn max for reclustering algorithms

#---YOUR DATASET TO READ IN
files   <- list.files()
files   <- files[grepl(".gml" ,files,fixed=T)]

studies <- unlist(strsplit(files,".gml"))

#--PPI networks
gg <- read.graph(files[1],format="gml")
gg <- simplify(gg,remove.multiple=T)

cat("Network: ", studies[1], "N: ", length(V(gg)), " E: ", length(E(gg)),"\n")

IDS <- V(gg)$name;
ids <- V(gg)$name;

#---subsampling scheme
if( type == 1 ){

    nr  <- ceiling( length(E(gg))*(mask/100) )
    ggM <- delete_edges(gg,sample(E(gg),nr))

    cat("Masking ", nr, " Edges\n")
    cat("Network: ", studies[1], "N: ", length(V(ggM)), " E: ", length(E(ggM)),"\n")

}

if( type == 2 ){

    nr  <- ceiling( length(V(gg))*(mask/100) )
    ggM <- delete_vertices(gg,sample(V(gg),nr))

    cat("Masking ", nr, " Vertices\n")
    cat("Network: ", studies[1], "N: ", length(V(ggM)), " E: ", length(E(ggM)),"\n")

}


#---Find Largest CC
ggLCC <- findLCC(ggM)
cat("LCC for: ", studies[1], "N: ", length(V(ggLCC)), " E: ", length(E(ggLCC)),"\n")
#---


#---build consensus file
cc       <- matrix(-1, ncol=3, nrow=length(V(gg)))
cc[,1]   <- V(gg)$name
cc[,2]   <- ifelse(cc[,1] %in% V(ggLCC)$name,cc[,1],-1)


#---Calculate Cnmin
cat("\n")
cat(Cnmin, "% of network (ggLCC) size = ", length(V(ggLCC)))

if( Cnmin > 0 ){
  Cnmin = floor( (Cnmin*length(V(ggLCC)))/100 )
} else {
  Cnmin = 1;
}
cat(", Cnmin set to: ", Cnmin, "\n")
#---

#---Calculate Cnmax
cat("\n")
cat(Cnmax, "% of network (ggLCC) size = ", length(V(ggLCC)))

if( Cnmax > 0 ){
  Cnmax = floor( (Cnmax*length(V(ggLCC)))/100 )
} else {
#---default is 10% of network size
  Cnmax = floor( (10*length(V(ggLCC)))/100 )
}
cat(", Cnmax set to: ", Cnmax, "\n")
#---


#---igraph's clustering algs.
alg     <- vector(length=16)
alg[1]  <- "louvain"
alg[2]  <- "infomap"
alg[3]  <- "fc"
alg[4]  <- "lec"
alg[5]  <- "sgG1"
alg[6]  <- "sgG2"
alg[7]  <- "sgG5"
alg[8]  <- "wt"
alg[9]  <- "spectral"
alg[10] <- "louvain2"
alg[11] <- "infomap2"
alg[12] <- "fc2"
alg[13] <- "lec2"
alg[14] <- "sgG12"
alg[15] <- "wt2"
alg[16] <- "Spectral2"

cat("running ", alg[aa], "...")

#lourvain
if( aa == 1 ) {

 louvain <- cluster_louvain(ggLCC)
 cc[,3]   <- ifelse(cc[,2] %in% louvain$names,louvain$membership,-1)

}

#infomap
if( aa == 2 ){

    infomap <- cluster_infomap(ggLCC)
    cc[,3]  <- ifelse(cc[,2] %in% infomap$names,infomap$membership,-1)

}


if( aa == 3 ){

    fc  <- fastgreedy.community(ggLCC)
    cc[,3]  <- ifelse(cc[,2] %in% fc$names,fc$membership,-1)

}


if( aa == 4 ){

    lec     <- leading.eigenvector.community(ggLCC)
    cc[,3]  <- ifelse(cc[,2] %in% lec$names,lec$membership,-1)
}


if( aa == 5 ){

    sg  <- spinglass.community(ggLCC, spins=as.numeric(500),gamma=1)
    cc[,3]  <- ifelse(cc[,2] %in% sg$names,sg$membership,-1)


}

if( aa == 6 ){

    sg  <- spinglass.community(ggLCC, spins=as.numeric(500),gamma=2)
    cc[,3]  <- ifelse(cc[,2] %in% sg$names,sg$membership,-1)


}


if( aa == 7 ){

    sg  <- spinglass.community(ggLCC, spins=as.numeric(500),gamma=5)
    cc[,3]  <- ifelse(cc[,2] %in% sg$names,sg$membership,-1)


}


if( aa == 8 ){

    wt  <- walktrap.community(ggLCC)
    cc[,3]  <- ifelse(cc[,2] %in% wt$names,wt$membership,-1)

}


if( aa == 9 ){

    el <- as.data.frame(get.edgelist(ggLCC,names=T))

    spec   <- CDMSuite::spectral(DF=el, CnMIN=Cnmin);
    cc[,3] <- ifelse(cc[,2] %in% spec$ID,spec$K,-1)

}


if( aa == 10 ){

    #---STEP 1
    #cluster with louvain, and recluster
    louvain <- cluster_louvain(ggLCC)

    ggLCC = igraph::set.vertex.attribute(ggLCC,"louvain",V(ggLCC),louvain$membership)

    oo = recluster( ggLCC, "louvain", Cnmax )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}

if( aa == 11 ){

    #---STEP 1
    #cluster with infomap, and recluster
    infomap <- cluster_infomap(ggLCC)

    ggLCC = igraph::set.vertex.attribute(ggLCC,"infomap",V(ggLCC),infomap$membership)

    oo = recluster( ggLCC, "infomap", Cnmax )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}


if( aa == 12 ){

    #---STEP 1
    #cluster with fc, and recluster
    fc <- fastgreedy.community(ggLCC)

    ggLCC = igraph::set.vertex.attribute(ggLCC,"fc",V(ggLCC),fc$membership)

    oo = recluster( ggLCC, "fc", Cnmax )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}

if( aa == 13 ){

    #---STEP 1
    #cluster with lec, and recluster
    lec     <- leading.eigenvector.community(ggLCC)

    ggLCC = igraph::set.vertex.attribute(ggLCC,"lec",V(ggLCC),lec$membership)

    oo = recluster( ggLCC, "lec", Cnmax )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}

if( aa == 14 ){

    #---STEP 1
    #cluster with sg, and recluster
    sg  <- spinglass.community(ggLCC, spins=as.numeric(500),gamma=1)

    ggLCC = igraph::set.vertex.attribute(ggLCC,"sgG1",V(ggLCC),sg$membership)

    oo = recluster( ggLCC, "sgG1", Cnmax )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}

if( aa == 15 ){

    #---STEP 1
    #cluster with sg, and recluster
    wt  <- walktrap.community(ggLCC)

    ggLCC = igraph::set.vertex.attribute(ggLCC,"wt",V(ggLCC),wt$membership)

    oo = recluster( ggLCC, "wt", Cnmax )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}

if( aa == 16 ){

    #---STEP 1
    #cluster with Spectral, and recluster
    el   <- as.data.frame(get.edgelist(ggLCC,names=T))
    spec <- CDMSuite::spectral(DF=el, CnMIN=Cnmin);

    mem   = spec$K[match(V(ggLCC)$name,spec$ID)]

    ggLCC = igraph::set.vertex.attribute(ggLCC,"Spectral",V(ggLCC), mem)

    oo = recluster( ggLCC, "Spectral", Cnmax, CnMIN=Cnmin )

    if( !is.null(oo) ){
        cc[,3]   <- ifelse(cc[,2] %in% oo[,1],oo[,4],-1)
    }

}


cat("done. \n")

cc <- as.data.frame(cc)
outfile <- file("consensusout.txt","w")
cat("#consensus",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
