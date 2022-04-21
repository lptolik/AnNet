##
library(igraph)
library(topGO)
library(org.Hs.eg.db)

## The GO Ontology to run the enrichment
TAX       = "9606"
Onto      = "CC"
TAXmap    = org.Hs.eg.db
TAXmapStr = "org.Hs.eg.db"    
TAXtitle  = "Human"

## Get Gene Symbols 
geneID2TERM <- revmap(topGO::annFUN.org(feasibleGenes=NULL,whichOnto=Onto, mapping = TAXmapStr, ID = "Entrez"))
TERM2geneID <- revmap(geneID2TERM)

geneIDs=names(geneID2TERM)
##---


##---Get Gene Entrez IDs from graph
#--- load corresponding graph 
gg = igraph::read.graph(sprintf("PPI_Presynaptic.gml",format="gml")
myList <- as.vector(V(gg)$name)
myList <- unique(myList)
myList <- na.omit(myList)
##---

## Find graph entrez IDs in GO ontology annotation set
geneList <- factor(as.integer(geneIDs %in% myList))
names(geneList) <- geneIDs

## setup the topGO enrichment object 
GOdata <- new("topGOdata", ontology=Onto, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2TERM)

##--- Run enrichment using topGO
resultFisher        <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher.elim   <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
        
topnodes = resultFisher@geneData[[4]]

if( resultFisher.elim@geneData[[4]] < topnodes ){
    topnodes = resultFisher.elim@geneData[[4]]
}

## The GO object after running the enrichment analysis
allRes <- GenTable(GOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topnodes )
##---

## Add Entrez IDs to output GO object
hits=lapply(allRes$GO.ID,function(x)
                          intersect(sigGenes(GOdata),genesInTerm(GOdata, c(x))[[1]]))
allRes$hits=hits
##---


## Create Flat file
## Note at this point you can pass a list of GO terms to include, or use instead,
## and increase/decrease the thresholds.
## 
thresMIN <- 10
thresMAX <- 1000

NGOterms = dim(allRes)[1]

flatfile <- data.frame()

for( i in 1:NGOterms ){

    ids  = allRes$hits[i][[1]]
    Nids = length(ids)
    if( Nids >= thresMIN && Nids <= thresMAX ){

        if( i == 1 ){
            flatfile = cbind(rep(allRes$GO.ID[i][[1]],Nids),
                             rep(allRes$Term[i][[1]], Nids),
                             ids)
        } else {
            tmp = cbind(rep(allRes$GO.ID[i][[1]],Nids),
                        rep(allRes$Term[i][[1]], Nids),
                        ids)
            flatfile = rbind(flatfile, tmp)
        }        

    }
    
}


write.table(flatfile,
            sprintf("flatfile.go.%s.csv",Onto),
            sep="\t",
            col.names=F,
            row.names=F,
            quote=F)
##---
