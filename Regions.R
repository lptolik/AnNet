#---Define Quadrant (or Regions) each gene lies in

#---Set default values for parameters
Nreg     <- NULL      #---Number of Regions in B V. SL plot in Bridgeness.R
Nele     <- NULL      #---Number of x,y elements needed for Scheme
quad     <- NULL
quadName <- NULL
Xmin     <- NULL
Xmax     <- NULL
Ymin     <- NULL
Ymax     <- NULL
#Xoff     <- NULL
MainDivSize <- 0.8

cat("***********\n")
cat("Loading Bridging Scheme: ", Scheme, " defined in 'Regions.R' \n")
cat("***********\n")

#Basic scheme 
if( Scheme == 1 ){

    Nreg <- 4
    Nele <- 4
    
    
    quad     <- vector(length=Nele)
    quadName <- vector(length=Nele)

    Xmin <- vector(length=Nele)
    Xmax <- vector(length=Nele)
    Ymin <- vector(length=Nele)
    Ymax <- vector(length=Nele)

    #quadrant each gene lies in
    Xmin[1] = 0.0; Xmax[1] = 0.5; Ymin[1] = 0.5; Ymax[1] = 1.0; quad[1] = 1; quadName[1] = "ul";
    Xmin[2] = 0.5; Xmax[2] = 1.0; Ymin[2] = 0.5; Ymax[2] = 1.0; quad[2] = 2; quadName[2] = "ur";
    Xmin[3] = 0.0; Xmax[3] = 0.5; Ymin[3] = 0.1; Ymax[3] = 0.5; quad[3] = 3; quadName[3] = "ll";
    Xmin[4] = 0.5; Xmax[4] = 1.0; Ymin[4] = 0.0; Ymax[4] = 0.5; quad[4] = 4; quadName[4] = "lr";

    Gseg1 <- geom_segment(aes(x=Xmin[1],xend=Xmax[1],y=0.1,yend=0.1),colour="grey40",size=MainDivSize,linetype=2,show.legend=F)
    Gseg2 <- geom_segment(aes(x=Xmin[1],xend=Xmin[1],y=Ymin[1],yend=Ymin[1]),colour="grey40",size=0,linetype=2,show.legend=F)
    Gseg3 <- geom_segment(aes(x=Xmin[1],xend=Xmin[1],y=Ymin[1],yend=Ymax[1]),colour="grey40",size=0,linetype=2,show.legend=F)
    

    MainDivSize = 1;    

    
    GAnno1 <- annotate("text",x=0.45, y=0.95,colour="limegreen",size=rel(15),label=quad[1])
    GAnno2 <- annotate("text",x=0.95, y=0.95,colour="limegreen",size=rel(15),label=quad[2])
    GAnno3 <- annotate("text",x=0.45, y=0.15,colour="limegreen",size=rel(15),label=quad[3])
    GAnno4 <- annotate("text",x=0.95, y=0.05,colour="limegreen",size=rel(15),label=quad[4])
    GAnno5 <- annotate("text",x=0.0,  y=0.0,colour="limegreen",size=0,label="")
}

#Deferential non-bridging, primary and secondray, bridging proteins using pathway annotation. 
if( Scheme == 2 ){

    Nreg <- 5
    Nele <- 6

    quad     <- vector(length=Nele)
    quadName <- vector(length=Nele)

    Xmin <- vector(length=Nele)
    Xmax <- vector(length=Nele)
    Ymin <- vector(length=Nele)
    Ymax <- vector(length=Nele)

    #quadrant each gene lies in
    Xmin[1] = 0;     Xmax[1] = 0.5;   Ymin[1] = 0;     Ymax[1] = 0.1;   quad[1] = 1; quadName[1] = "Non-Bridging";
    Xmin[2] = 0;     Xmax[2] = 0.125; Ymin[2] = 0.1;   Ymax[2] = 0.375; quad[2] = 2; quadName[2] = "Bridging-Synaptic-PW";
    Xmin[3] = 0.125; Xmax[3] = 1.0;   Ymin[3] = 0.1;   Ymax[3] = 0.5;   quad[3] = 3; quadName[3] = "Secondary-Bridging-glob";
    Xmin[4] = 0;     Xmax[4] = 0.125; Ymin[4] = 0.375; Ymax[4] = 0.5;   quad[4] = 4; quadName[4] = "Bridging-Immune-PW";
    Xmin[5] = 0;     Xmax[5] = 0.125; Ymin[5] = 0.5;   Ymax[5] = 1.0;   quad[5] = 4; quadName[5] = "Bridging-Immune-PW";
    Xmin[6] = 0.125; Xmax[6] = 1.0;   Ymin[6] = 0.5;   Ymax[6] = 1.0;   quad[6] = 5; quadName[6] = "Primary-Bridging-glob";

    Gseg1 <- geom_segment(aes(x=Xmin[1],xend=Xmax[6],y=Ymax[1],yend=Ymax[1]),colour="grey40",size=1,linetype=2,show.legend=F)
    Gseg2 <- geom_segment(aes(x=Xmin[1],xend=Xmax[6],y=Ymin[4],yend=Ymin[4]),colour="grey40",size=1,linetype=2,show.legend=F)
    Gseg3 <- geom_segment(aes(x=Xmax[2],xend=Xmax[2],y=Ymin[2],yend=Ymax[6]),colour="grey40",size=1,linetype=2,show.legend=F)
    
    MainDivSize = 0;

    GAnno1 <- annotate("text",x=0.95,  y=0.05, colour="limegreen",size=rel(15),label=quad[1])
    GAnno2 <- annotate("text",x=0.095, y=0.3,  colour="limegreen",size=rel(14),label=quad[2])
    GAnno3 <- annotate("text",x=0.95,  y=0.3,  colour="limegreen",size=rel(15),label=quad[3])
    GAnno4 <- annotate("text",x=0.05,  y=0.95, colour="limegreen",size=rel(15),label=quad[4])
    GAnno5 <- annotate("text",x=0.95,  y=0.95, colour="limegreen",size=rel(15),label=quad[6])
    
}

#Entropy Scheme 
if( Scheme == 3 ){

    Nreg <- 4
    Nele <- 4

    quad     <- vector(length=Nele)
    quadName <- vector(length=Nele)

    Xmin <- vector(length=Nele)
    Xmax <- vector(length=Nele)
    Ymin <- vector(length=Nele)
    Ymax <- vector(length=Nele)

    #quadrant each gene lies in
    Xmin[1] = 0.0;   Xmax[1] = IntCon; Ymin[1] = 0.5; Ymax[1] = 1.0; quad[1] = 1; quadName[1] = "ul";
    Xmin[2] = IntCon; Xmax[2] = 1.0;   Ymin[2] = 0.5; Ymax[2] = 1.0; quad[2] = 2; quadName[2] = "ur";
    Xmin[3] = 0.0;   Xmax[3] = IntCon; Ymin[3] = 0.1; Ymax[3] = 0.5; quad[3] = 3; quadName[3] = "ll";
    Xmin[4] = IntCon; Xmax[4] = 1.0;   Ymin[4] = 0.0; Ymax[4] = 0.5; quad[4] = 4; quadName[4] = "lr";

    Gseg1 <- geom_segment(aes(x=Xmin[1],xend=Xmax[1],y=0.1,yend=0.1),colour="grey40",size=MainDivSize,linetype=2,show.legend=F)
    Gseg2 <- geom_segment(aes(x=Xmin[1],xend=Xmin[1],y=Ymin[1],yend=Ymin[1]),colour="grey40",size=0,linetype=2,show.legend=F)
    Gseg3 <- geom_segment(aes(x=Xmin[1],xend=Xmin[1],y=Ymin[1],yend=Ymax[1]),colour="grey40",size=0,linetype=2,show.legend=F)    

    MainDivSize = 1;    

    
    GAnno1 <- annotate("text",x=(IntCon-0.05), y=0.95,colour="limegreen",size=rel(15),label=quad[1])
    GAnno2 <- annotate("text",x=0.95,          y=0.95,colour="limegreen",size=rel(15),label=quad[2])
    GAnno3 <- annotate("text",x=(IntCon-0.05), y=0.15,colour="limegreen",size=rel(15),label=quad[3])
    GAnno4 <- annotate("text",x=0.95,          y=0.05,colour="limegreen",size=rel(15),label=quad[4])
    GAnno5 <- annotate("text",x=0.0,           y=0.0, colour="limegreen",size=0,label="")
}


#
fillRegions <- function( DF, X, Y, indx ){


    if( Scheme == 1 ){

        #quadrant each gene lies in
        DF[which( X >= Xmin[1] & X <  Xmax[1] & Y >= Ymin[1] & Y <= Ymax[1]),indx] <- quad[1] #ul
        DF[which( X >= Xmin[2] & X <= Xmax[2] & Y >= Ymin[2] & Y <= Ymax[2]),indx] <- quad[2] #ur
        DF[which( X >= Xmin[3] & X <  Xmax[3] & Y >= Ymin[3] & Y <  Ymax[3]),indx] <- quad[3] #ll
        DF[which( X >= Xmin[4] & X <= Xmax[4] & Y >= Ymin[4] & Y <  Ymax[4]),indx] <- quad[4] #lr
        
    }
    
    
    if( Scheme == 2 ){        

        #---Quadrant (or Region) each gene lies in
        DF[which( X >= Xmin[1] & X <  Xmax[1] & Y >= Ymin[1] & Y <  Ymax[1]),indx] <- quad[1] #non
        DF[which( X >= Xmin[2] & X <  Xmax[2] & Y >= Ymin[2] & Y <  Ymax[2]),indx] <- quad[2] #sec-bridge-loc
        DF[which( X >= Xmin[3] & X <  Xmax[3] & Y >= Ymin[3] & Y <  Ymax[3]),indx] <- quad[3] #sec-bridge-loc-glob
        DF[which( X >= Xmin[4] & X <  Xmax[4] & Y >= Ymin[4] & Y <  Ymax[4]),indx] <- quad[4] #prim-bridge-loc-2
        DF[which( X >= Xmin[5] & X <  Xmax[5] & Y >= Ymin[5] & Y <= Ymax[5]),indx] <- quad[5] #prim-bridge-loc-1        
        DF[which( X >= Xmin[6] & X <= Xmax[6] & Y >= Ymin[6] & Y <= Ymax[6]),indx] <- quad[6] #prim-bridge-glob

    }
        
    if( Scheme == 3 ){

        #quadrant each gene lies in
        DF[which( X >= Xmin[1] & X <  Xmax[1] & Y >= Ymin[1] & Y <= Ymax[1]),indx] <- quad[1] #ul
        DF[which( X >= Xmin[2] & X <= Xmax[2] & Y >= Ymin[2] & Y <= Ymax[2]),indx] <- quad[2] #ur
        DF[which( X >= Xmin[3] & X <  Xmax[3] & Y >= Ymin[3] & Y <  Ymax[3]),indx] <- quad[3] #ll
        DF[which( X >= Xmin[4] & X <= Xmax[4] & Y >= Ymin[4] & Y <  Ymax[4]),indx] <- quad[4] #lr
        
    }
    
    
    return(DF)
    
}
