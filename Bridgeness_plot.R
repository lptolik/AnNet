source('../setUp.R');

scale <- function(x, VALUE=NULL){

    x = as.numeric(as.vector(x))

    xmin <- min(x,na.rm=T)
    xmax <- max(x,na.rm=T)
   
    if( is.null(VALUE) ){

        x  <- x-xmin
        x  <- ifelse(!is.na(x), x/(xmax-xmin), NA) 

        return(x)
    }

    value = as.numeric(as.vector(value)[1])
    value = value-xmin
    value = ifelse(!is.na(value), value/(xmax-xmin), NA) 
    return(value)
}

##--- Directories needed
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]

    
##--- declare clustering algorithms used 
##--- Clustering algorithms used 
alg <- "Spectral"

##--- load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

##--- No: nodes in graph
N  <- length(V(gg))

##--- read Semi-local Centrality 'Cl' from graph 'gg'.
X    <- igraph::get.vertex.attribute(gg,"Cl",V(gg))
X    <- scale(X)
Xlab <- "Local Centrality (Cl)" 
Ylab <- "Bridgeness (B)"

##---And set Scheme in 'Regions.R' for a Bridgeness V. Semi-local Centrality 'Cl' Plot
Scheme      <- 1
MainDivSize <- 0.8
xmin        <- 0
xmax        <- 1
ymin        <- 0
ymax        <- 1
source("Regions.R")
##----

    
##---Check or create output plot dir
plotDIR <- sprintf("%s/PLOTS/",subDIR[S])

if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}

##--- set random no: seed and for geom_text_repel
set.seed(42)

#Force1    <- vector(length=3)
#Force1[1] <- 10
#Force1[2] <- 5
#Force1[3] <- 10

#Force3    <- vector(length=3)
#Force3[1] <- 10
#Force3[2] <- 1
#Force3[3] <- 10
##----    

##--- Read-indx gene lists
##    Note: these gene list should be added as a vertex attribute and
##          read-in from the graph 'gg'

##--- Read-in Core PSD95 interactor genes
CORE  <- read.table("CorePSD95Complex.csv",sep="\t",header=F)[[1]]
indx  <- match(V(gg)$name,CORE)
group <- ifelse( is.na(indx), 0,1)
    
##--- Read-in my 'VIP' gene list
VIPs   <- read.table("BrProteins_2020.csv",sep="\t",header=F)[[1]]
indx   <- match(V(gg)$name,VIPs)
group2 <- ifelse( is.na(indx), 0,1)
##---
    
##---
## Build Bridging (B) V. Semi-local Centrality (Cl) plot 
##---

cat("Plotting Br V. ", Xlab ," for: ", alg, "\n")

alg.br  <- sprintf("BRIDGE_%s",alg)   
bridge  <- igraph::get.vertex.attribute(gg,alg.br,V(gg))
Y       <- as.numeric(as.vector(bridge))    


##--- container to hold Quadrant (or Regions) each gene lies in
##    given each genes Bridgesness and centrality value
cons <- matrix(0,ncol=3,nrow=N)
colnames(cons) <- c("EntezID","GeneName","reg")
cons[,1] <- igraph::get.vertex.attribute(gg,"name",V(gg))     ## Entrez ID
cons[,2] <- igraph::get.vertex.attribute(gg,"GeneName",V(gg)) ## Gene Symbol
cons[,3] <- rep(0,N)

##--- Fill quadrant (or Region) each gene lies in
##    'fillRegions' defined in 'Regions.R' 
cons[,3] <- rep(0,length(cons[,1])) 
cons     <- fillRegions( cons, X, Y, 3 );      

##--- First-Level GeneNames
##    record each gene name in each region
genes1 <- ifelse( cons[,3]==quad[1],cons[,2],"")
genes2 <- ifelse( cons[,3]==quad[2],cons[,2],"")
genes3 <- ifelse( cons[,3]==quad[3],cons[,2],"")
genes4 <- ifelse( cons[,3]==quad[4],cons[,2],"")
genes5 <- ifelse( cons[,3]==quad[5],cons[,2],"")
genes6 <- ifelse( cons[,3]==quad[6],cons[,2],"")

##--- Second-Level GeneNames
##    Filter First-Level GeneNames using our VIP gene list
##    Note: We'll keep the First-Level GeneNames for UR region
VIPsGN <- cons[match(VIPs,cons[,1]),2]   
      
genes1 <- ifelse(!is.na(match(genes1,VIPsGN)),genes1,"")
##genes2 <- ifelse(!is.na(match(genes2,VIPsGN)),genes2,"")
genes3 <- ifelse(!is.na(match(genes3,VIPsGN)),genes3,"")
genes4 <- ifelse(!is.na(match(genes4,VIPsGN)),genes4,"")
genes5 <- ifelse(!is.na(match(genes5,VIPsGN)),genes5,"")
genes6 <- ifelse(!is.na(match(genes6,VIPsGN)),genes6,"")        

## GeneName to plot in:
## UL = Upper Left
## UR = Upper Right
## LL = Lower Left
## LR = Lower Right
GenesUL <- genes1
GenesUR <- genes2
GenesLL <- genes3
GenesLR <- genes4

#if( Scheme == 2 ){            
#    GenesUL <- genes4
#    GenesUR <- genes6
#    GenesLL <- genes2
#    GenesLR <- genes3
#}

##--- Build Br V Cl ggplot data.frame
df      = cbind(X,Y,group,group2)
df      = as.data.frame(df)
df_core = df[df[,3]==1,]
df_sp   = df[df[,4]==1 & df[,3] !=1,]
df_bs   = df[df[,3] != 1 & df[,4] != 1,]       

##--- baseColor of genes
baseColor="royalblue2"

##--- SPcolor, colour highlighting any 'specical' genes
SPColor="royalblue2"

##--- PSDColor, colour of core PSD95 interactor genes
PSDColor="magenta"
        
##--- ggplot object for Br V Cl plot
gplot <- ggplot(df,aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)) ))+                          
    ##base
    geom_point(data=df_bs,
               aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)), alpha=(X*Y)),colour=baseColor,shape=16,show.legend=F)+
    ##Special
    geom_point(data=df_sp,
               aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y))),colour=SPColor,alpha=1,shape=16,show.legend=F)+
    ##core PSD95
    geom_point(data=df_core,
               aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y))),colour=PSDColor,alpha=1,shape=16,show.legend=F)+
    ##ul
    geom_label_repel(aes(label=as.vector(GenesUL)),
                     fontface='bold',color='black',fill='white',box.padding=0.1,point.padding=NA,label.padding=0.15,segment.color='black',
                     xlim=c(0.0,0.5), ylim=c(0.5,1.0),force=1,size=rel(3.8),show.legend=F)+
    ##ur
    geom_text_repel(aes(label=as.vector(GenesUR)),
                    fontface='bold',color='black',segment.color='grey50',point.padding=NA,
                    xlim=c(0.45,1.0), ylim=c(0.45,1.0),force=0.5,size=rel(3.8),show.legend=F)+  
    ##ll (i.e. prim-loc)
    geom_label_repel(aes(label=as.vector(GenesLL)),
                     fontface="bold",color="black",fill='white',box.padding=0.0,point.padding=NA,label.padding=0.1,segment.color='black',
                     xlim=c(0.0,0.5), ylim=c(0.1,0.5),force=0.5,size=rel(3.0),show.legend=F)+
    ##lr
    geom_text_repel(aes(label=as.vector(GenesLR)),
                    color="black",point.padding=NA,
                    xlim=c(0.5,1.0), ylim=c(0.0,0.5),force=5,size=rel(3.8),show.legend=F)+ 
    
    labs(x=Xlab,y=Ylab,title=sprintf("%s",alg))+
    scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax))+
    theme(            
        axis.title.x=element_text(face="bold",size=rel(2.5)),
        axis.title.y=element_text(face="bold",size=rel(2.5)),
        legend.title=element_text(face="bold",size=rel(1.5)),
        legend.text=element_text(face="bold",size=rel(1.5)),
        legend.key=element_blank())+
    theme(panel.grid.major = element_line(colour="grey40",size=0.2),
          panel.grid.minor = element_line(colour="grey40",size=0.1),
          panel.background = element_rect(fill="white"),
          panel.border = element_rect(linetype="solid",fill=NA))+
    
    geom_vline(xintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
    geom_hline(yintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+

    ## GsegX and GAnnoX defined in 'Regions.R' 
    Gseg1 +
    Gseg2 +
    Gseg3 +
    
    GAnno1 +
    GAnno2 +
    GAnno3 +
    GAnno4 +
    GAnno5

png(sprintf("%s/%s_%s_BrV%s.png",plotDIR,subDIR[S],alg[a],Xmeas),width=WIDTH,height=HEIGHT,units="px")
print(gplot)
dev.off()
    
