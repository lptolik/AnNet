#' Calculate annotation enrichment for clusters in the graph
#'
#' @param g graph to get annotation from
#' @param alg cluster algorithm and membership attribute name
#' @param name annotation attribute name
#' @param col list separation character in attribute, by
#' default is \code{;}
#' @param vid attribute to be used as a vertex ID
#' @param alpha probability threshold
#'
#' @return A table with ORA results.
#' Each row corresponds to a tested annotation in particular cluster.
#' The columns are the following:
#' \itemize{
#'   \item pathway – name of the pathway as in 'names(pathway)';
#'   \item pval – an enrichment p-value from hypergeometric test;
#'   \item padj – a BH-adjusted p-value;
#'   \item overlap – size of the overlap;
#'   \item size – size of the gene set;
#'   \item leadingEdge – vector with overlapping genes.
#'   \item cl – cluster ID
#'   }
#' @export
#' @importFrom fgsea fora
clusterORA<-function(g,alg,name,vid='name',alpha=0.1,col=COLLAPSE){
  anL<-getAnnotationVertexList(g,'GOCCID')
  cl<-make_clusters(g,as.numeric(get.vertex.attribute(g,alg)))
  forafun<-function(.i){
    res<-as.data.frame(fora(anL,
              get.vertex.attribute(g,vid)[which(membership(cl)==.i)],
              universe = as.character(get.vertex.attribute(g,vid))));
    res$cl<-.i
    return(res)
  }
  resL<-lapply(1:length(cl),forafun)
  res<-do.call(rbind,resL)
  res<-res[res$padj<alpha,]
  res$overlapGenes<-sapply(res$overlapGenes,paste,collapse = ', ')
  return(res)
}
