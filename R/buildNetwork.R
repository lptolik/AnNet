#---Find Largest CC
#' Find Largest CC
#'
#' @param GG graph to analyze
#'
#' @return igraph representation largest CC
#' @export
#' @import igraph
#'
#' @examples
findLCC <- function(GG){

  dec <- decompose.graph(GG)
  d=1
  CC=length(V(dec[[1]]))
  for( i in 1:length(dec) ){
    if(length(V(dec[[i]])) > CC){
      d=i
      CC=length(V(dec[[i]]))
    }
  }

  GG  <- decompose.graph(GG)[[d]]
  return(GG)

}

#' Title
#'
#' @param eatt
#' @param TERMS
#'
#' @return
#' @export
#'
#' @examples
findTERM <- function(eatt, TERMS){

  eatt  <- as.vector(eatt)
  found <- rep(FALSE, length(eatt))

  T = length(TERMS)
  for( t in 1:T ){
    temp  <- grepl(eatt[t], eatt)
    found <- as.logical(found) | as.logical(temp)
  }

  return(found)

}


#' Title
#'
#' @param TERMS
#'
#' @return
#' @export
#'
#' @examples
filterPMIDs <- function(TERMS=NULL){

  filterIDs <- c()

  if( !is.null(TERMS ) && length(TERMS) != 0 ){

    pmids <- read.delim("/afs/inf.ed.ac.uk/user/c/cmclean5/ownCloud/Synaptic_proteome/anaysis_17_05_2019/mined_PPIs/pmid_keywords.csv",sep="\t",header=T)

    indX <- list()

    for( i in 1:length(TERMS) ){
      N = length(names(indX))
      indX[[N+1]] <- grepl(TERMS[i], pmids[,4])
      names(indX)[N+1] <- sprintf("%s_title",TERMS[i])
      N = length(names(indX))
      indX[[N+1]] <- grepl(TERMS[i], pmids[,5])
      names(indX)[N+1] <- sprintf("%s_keywords",TERMS[i])
    }

    exc <- indX[[1]]
    for( i in 2:length(names(indX)) ){
      exc <- exc | indX[[i]]
    }

    filterIDs <- as.vector(pmids[exc,1])

  }

  return(filterIDs)

}

#--- add attributes to igraph edges from it raw file
#' Title
#'
#' @param GG
#' @param gg
#'
#' @return
#' @export
#'
#' @examples
addEdgeAtts <- function(GG, gg){

  ATTS = names(edge.attributes(GG))

  if( !is.null(ATTS) ){

    ed = get.edgelist(gg)
    M  = length(E(gg))
    ED = get.edgelist(GG)

    VALUES = list()

    for( a in 1:length(ATTS) ){
      VALUES[[a]] = get.edge.attribute(GG,ATTS[a],E(GG))
      names(VALUES)[a] = ATTS[a]
    }

    cat("\n")
    cat("scanning edges...")
    RES    = matrix("",nrow=M, ncol=length(ATTS))

    for( e in 1:M ){

      indx = (ed[e,1] == ED[,1] & ed[e,2] == ED[,2]) | (ed[e,1] == ED[,2] & ed[e,2] == ED[,1])

      for( a in 1:length(ATTS) ){

        res = VALUES[[a]][indx]

        if( res != "" ){
          res <- unique(res)
          if( length(res) == 1 ){
            RES[e,a] <- res
          } else {
            RES[e,a] <- paste(as.character(res),collapse=';')
          }
        }
      }
    }

    cat("done.\n")

    for( a in 1:length(ATTS) ){
      gg <- set.edge.attribute(gg,ATTS[a],E(gg),as.character(RES[,a]))
    }

  }

  return(gg)

}

#' Build network from data.table
#'
#' @param ff network structure data.frame with first two columns defining the
#' network edge nodes
#' @param kw pmid keyword annotation data.frame. If `NA`
#' no annotation will be added
#'
#' @return igraph object of the largest connected component
#' @export
#' @import igraph
#'
#' @examples
#' f<-data.frame(A=c('A','A','B'),B=c('B','C','C'))
#' gg<-buildNetwork(f)
buildNetwork<-function(ff,kw=NA){
  #--- build raw graph
  GG <- graph.data.frame(ff[,1:2],directed=F)
  if( !is.na(kw) ){
    GG = set.edge.attribute(GG,"METHOD",E(GG), as.character(ff[,3]))
    GG = set.edge.attribute(GG,"TYPE",E(GG), as.character(ff[,7]))

    PMIDS = ifelse(!grepl("unassigned",ff[,4]), sprintf("PMID:%s",ff[,4]), ff[,4])
    GG = set.edge.attribute(GG,"PUBMED",E(GG), PMIDS)

    YEARS = kw[match(gsub("PMID:","",E(GG)$PUBMED),kw[,1]),3]
    YEARS = ifelse(is.na(YEARS),"na",YEARS)
    GG = set.edge.attribute(GG,"YEAR",E(GG), YEARS)
    #---

  }
  #--- build igraph, removing multiple edges and loops
  gg <- simplify(GG,remove.multiple=T,remove.loops=T)
  #---Find Largest CC
  gg  <- findLCC(gg)
  gg <- addEdgeAtts(GG,gg)
}

#' Title
#'
#' @return
#' @export
#' @import synaptome.db
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
buildFromSynaptomeByEntrez<-function(entrez){
  t<-findGenesByEntrez(entrez)
  gg<-buildFromSynaptomeGeneTable(t)
  return(gg)
}

#' Title
#'
#' @param t
#'
#' @return
#' @export
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t)
buildFromSynaptomeGeneTable<-function(t){
  p<-getPPIbyIDs(t$GeneID,type = 'limited')
  aidx<-match(p$A,t$GeneID)
  bidx<-match(p$B,t$GeneID)
  gg<-buildNetwork(data.frame(A=t$HumanEntrez[aidx],B=t$HumanEntrez[bidx]))
  return(gg)
}

#' Calculate sparsness of the graph.
#'
#' @param gg graph to evaluate
#'
#' @return sparsness value
#' @export
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t)
#' calcSparsness(gg)
calcSparsness<-function(gg){
  N<-igraph::vcount(gg)
  E<-igraph::ecount(gg)
  sp<-2.0*E/(N*(N-1))
  return(sp)
}
