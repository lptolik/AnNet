#run 'mergeRESULTS.R' first

source('../setUp.R')
require(cowplot)

makeBold <- function(src, bolder) {
    require(purrr)                                            # make sure we've got purrr library
    if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
    src_levels <- levels(src)                                 # retrieve the levels in their order
    temp <- bolder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
    if (all(temp)) {                                         # if so
        b_pos <- purrr::map_int(bolder, ~which(.==src_levels)) # then find out where they are
        b_vec <- rep("plain", length(src_levels))               # make'm all plain first
        b_vec[b_pos] <- "bold"                                  # make our targets bold
        b_vec                                                   # return the new vector
    } else {
        stop("All elements of 'bolder' must be in src")
    }
}


FILES    <- vector(length=1)
FILES[1] <- "Disease_overlap_sig"

#--- load results file from calRANDOMDiseasePairs.R
ff <- read.table(sprintf("PPI_Presynaptic/%s.csv",subDIR[S],FILES[1]),sep="\t", header=T, stringsAsFactors=F, quote="", check.names=F)

## index of parameters of interest in result file
sabIDX <- which(grepl("sAB",colnames(ff),fixed=T)==TRUE)

zsIDX  <- which(grepl("zScore",colnames(ff),fixed=T)==TRUE)

BonIDX <- which(grepl("Bonferroni",colnames(ff),fixed=T)==TRUE)

pvIDX  <- which(grepl("pvalue",colnames(ff),fixed=T)==TRUE)

adIDX  <- which(grepl("p.adjusted",colnames(ff),fixed=T)==TRUE)
    
qvIDX  <- which(grepl("q-value",colnames(ff),fixed=T)==TRUE)
    
DDpairs <- sprintf("%s-%s",ff[,3],ff[,7])

## save parameters of interest
oo     <- matrix(0, ncol=7, nrow=length(DDpairs) )
colnames(oo) <- c("pairs","sAB","zscore","pvalue","p.adjusted","qvalue","label")
oo[,1] <- as.character(DDpairs)
oo[,2] <- as.character(ff[,sabIDX[1]])
oo[,3] <- as.character(ff[,zsIDX[1]])
oo[,4] <- as.character(ff[,pvIDX[1]])
oo[,5] <- as.character(ff[,adIDX[1]])
oo[,6] <- as.character(ff[,qvIDX[1]])
oo[,7] <- rep("Pre",length(DDpairs))

## remove disease "AUT" from plot
oo <- oo[oo[,2] != 0 & oo[,3] != 0,]
oo <- oo[!grepl("AUT",oo[,1]),]

## convert data into data.frame
df <- as.data.frame(oo)
df <- df[df$zscore < 0,] ## select only overlapping disease-pairs for plotting

## reorder parameters according to the pvalue 
ReOrder = order(as.numeric(df$pvalue),decreasing=T)
df = df[ReOrder,]   
df$pairs      <- as.vector(factor(df$pairs))
df$sAB        <- as.vector(factor(df$sAB))
df$zscore     <- as.vector(factor(df$zscore))
df$pvalue     <- as.vector(factor(df$pvalue))
df$p.adjusted <- as.vector(factor(df$p.adjusted))
df$qvalue     <- as.vector(factor(df$qvalue))
df$pairs      <- factor(df$pairs, levels=df$pairs)

## disease pairs to highlight
target <- c("AD-HTN","HTN-PD","AD-PD")

## colours & shapes                                 
## Presynaptic ==> firebrick2
## PSP         ==> royalblue2
## PSPc        ==> green2   
colours <- c('firebrick2','royalblue2','green2')
shapes  <- c(16,17,15)
    
## x limits for qvalues
XMIN= 0
XMAX=-log10(min(as.numeric(df$qvalue)))
if( XMAX < 50 ){ XMAX = 50 } 
    
##plot qvalues
 gplot <- ggplot(df,aes(y=as.factor(df$pairs),x=-log10(as.numeric(as.vector(df$qvalue))) ))+
     geom_point(alpha=I(0.8),colour=colours[1], shape=shapes[1], size=rel(3),
                show.legend=F)+
     scale_size_identity()+
     labs(y="Disease pairs",x="-log10(q-value)")+
     theme(legend.key=element_blank())+
     geom_vline(xintercept=-log10(0.05),colour="grey25",linetype="longdash",
                alpha=0.5, size=rel(1),show.legend=F)+
     theme_bw()+
     scale_x_continuous(limit = c(XMIN, XMAX))+
     theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"))+
     theme(axis.text.x = element_text(angle = 40, hjust = 1, face="bold",
                                      size=rel(1.5)),
              axis.text.y = element_text(face=makeBold(df$pairs,target),size=rel(1.5)),
              axis.title.x=element_text(face="bold",size=rel(2.0)),
              axis.title.y=element_text(face="bold",size=rel(2.0)),
              legend.title=element_text(face="bold",size=rel(1.7)),
              legend.text=element_text(face="bold",size=rel(1.7)),
              legend.position="top"
              )+
     guides(colour = guide_legend(override.aes = list(shape = shapes[1],
                                                      size=rel(7))),
               size   = FALSE,
               shape  = FALSE)+
     scale_color_manual("Compartment",breaks="Pre",values=colours[1])

    png("disease_pairs_overlaps.png",width=WIDTH,height=2*HEIGHT,units="px");
    print(gplot)
    dev.off()

