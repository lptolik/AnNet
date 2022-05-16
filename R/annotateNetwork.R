#TODO: convert all annotate_ functions to use either loopOverFiles, or
#      ontologies, or data.frames

#' Remove vertex property.
#'
#' @param GG igraph object
#' @param NAME name of the vertex property to remove
#'
#' @return
#' @export
#'
#' @examples
removeVertexTerm <- function(GG,NAME){

  if( !is.null(get.vertex.attribute(GG,NAME)) ){
    GG <- remove.vertex.attribute(GG,name=NAME)
  }

  if( !is.null(get.vertex.attribute(GG,gsub("_","",NAME))) ){
    GG <- remove.vertex.attribute(GG,name=gsub("_","",NAME))
  }

  return(GG)

}

COLLAPSE <- vector(length=2)
COLLAPSE[1] <- ";"
COLLAPSE[2] <- "&"
c=1


#' Annotate graph from list of files
#'
#' @param GG igraph object
#' @param FILES list of file path strings to read annotation from
#' @param NAME name of the vertex property
#' @param IDS vertex IDs
#' @param addIDS if TRUE NAME_ID property will be added
#'
#' @return
#' @export
#'
#' @examples
loopOverFiles <- function(GG, FILES, NAME, IDS, addIDS){

  for( f in 1:length(FILES) ){

    GG <- removeVertexTerm(GG, NAME[f])

    if( addIDS ){
      GG <- removeVertexTerm(GG, sprintf("%s_ID",NAME[f]))
    }

    if( file.exists(sprintf("./%s",FILES[f])) ){

      #--- Set Disease (geneRIF db) attributes in .gml graph
      set.vertex.attribute(GG, NAME[f],V(GG),"")

      annoF    <- read.table(FILES[f],sep="\t",skip=1,strip.white=T,quote="")
      annoFIDS <- as.character(annoF[,3])

      typeF   <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))

      oo     <- matrix("",ncol=3,nrow=length(IDS))
      oo[,1] <- IDS

      for( i in 1:length(IDS) ){

        ind1 = which(annoFIDS==IDS[i])

        Str1 <- "";
        Str2 <- "";

        if( length(ind1) != 0 ){

          if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
          else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

          if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
          else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }

          if( grepl(COLLAPSE[c], Str1) ){
            Str1 <- strsplit(Str1,COLLAPSE[c])[[1]]
            Str1 <- unique(Str1)
            if( length(Str1) > 1 ){
              Str1 <- paste(as.character(Str1),collapse=COLLAPSE[c])
            }
          }

          if( grepl(COLLAPSE[c], Str2) ){
            Str2 <- strsplit(Str2,COLLAPSE[c])[[1]]
            Str2 <- unique(Str2)
            if( length(Str2) > 1 ){
              Str2 <- paste(as.character(Str2),collapse=COLLAPSE[c])
            }
          }

          oo[i,2] <- Str1
          oo[i,3] <- Str2

        }

      }

      GG <- set.vertex.attribute(GG,NAME[f],V(GG),as.character(oo[,2]))

      if( addIDS ){
        GG <- set.vertex.attribute(GG,sprintf("%s_ID",NAME[f]),V(GG),as.character(oo[,3]))
      }


    }
  }


  return(GG)

}


#Add GeneNames
#' Title
#'
#' @param gg
#'
#' @return
#' @export
#' @import org.Hs.eg.db
#'
#' @examples
#' file <- system.file("extdata", "PPI_Presynaptic.gml", package = "AnNet")
#' gg <- igraph::read.graph(file,format="gml")
#' agg<-annotateGeneNames(gg)
annotateGeneNames<-function(gg){
  ids = V(gg)$name

  gn <- mapIds(org.Hs.eg.db,ids,column="SYMBOL",keytype="ENTREZID")

  gg <- removeVertexTerm(gg,"GeneName")

  set.vertex.attribute(gg,"GeneName",V(gg),"")
  V(gg)$GeneName = gn
  return(gg)
}

getDType<-function(){
  #---HDO Disease short names
  dtype  <- vector(length=12);
  dtype[1]   = "AD";
  dtype[2]   = "BD";
  dtype[3]   = "AUT";
  dtype[4]   = "SCH";
  dtype[5]   = "ASD";
  dtype[6]   = "Epi";
  dtype[7]   = "ID";
  dtype[8]   = "HTN";
  dtype[9]   = "HD";
  dtype[10]  = "PD";
  dtype[11]  = "FTD";
  dtype[12]  = "MS";
  #dtype[12]  = "DMH";
  #dtype[13]  = "CNSD";

  return(dtype)
}
#' Title
#'
#' @return vector of diseas IDs of interest
#' @export
#'
#' @examples
getDiseases<-function(){
  #---HDO ID DISEASES of INTEREST
  disn    <- vector(length=12);
  disn[1]  <- "DOID:10652"#Alzheimer's_disease"
  disn[2]  <- "DOID:3312"#bipolar_disorder"
  disn[3]  <- "DOID:12849"#autistic_disorder"
  disn[4]  <- "DOID:5419"#schizophrenia"
  disn[5]  <- "DOID:0060041"#autism_spectrum_disorder
  disn[6]  <- "DOID:1826"#epilepsy_syndrome
  disn[7]  <- "DOID:1059"
  disn[8]  <- "DOID:10763"
  disn[9]  <- "DOID:12858"
  disn[10] <- "DOID:14330"
  disn[11] <- "DOID:9255"
  disn[12] <- "DOID:2377"
  return(disn)
}
#Add topOnto_ovg
#' Title
#'
#' @param gg
#' @param dis
#'
#' @return
#' @export
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' afile<-system.file("extdata", "flatfile_human_gene2HDO.csv", package = "AnNet")
#' dis    <- read.table(afile,sep="\t",skip=1,header=F,strip.white=T,quote="")
#' #agg<-annotate_topOnto_ovg(gg,dis)
annotate_topOnto_ovg<-function(gg,dis){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"TopOnto_OVG")
  gg <- removeVertexTerm(gg,"TopOnto_OVG_HDO_ID")

  #--- Set Disease (geneRIF db) attributes in .gml graph
  set.vertex.attribute(gg,"TopOnto_OVG",V(gg),"")
  set.vertex.attribute(gg,"TopOnto_OVG_HDO_ID",V(gg),"")

  #par    <- read.table("flatfile_human_gene2HDO.parentTerm.csv",sep="\t",skip=1,strip.white=T,quote="")
  #dis    <- read.table("flatfile_human_gene2HDO.csv",sep="\t",skip=1,header=F,strip.white=T,quote="")

  #dis    <- rbind(dis,par)

  disIDS <- dis[,3]

  disn <-getDiseases()
  dtype<- getDType()

  for( i in 1:length(ids) ){

    ind1 = which(disIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){
#TDOD: refactor this code to work without disn
      disv <- as.vector(dis[ind1,1]);

      indx <- match(disv,disn)

      for( j in 1:length(disv) ){

        if( !is.na(indx[j]) ){

          if( Str1 == "" ) { Str1 <- as.character(dtype[indx[j]]) }
          else {
            Str1 <- paste(c(Str1,as.character(dtype[indx[j]])),collapse=COLLAPSE[c]) }

          if( Str2 == "" ) { Str2 <- as.character(disn[indx[j]]) }
          else {
            Str2 <- paste(c(Str2,as.character(disn[indx[j]])),collapse=COLLAPSE[c]) }
        }

      }
    }

    Str1 = paste(unique(strsplit(Str1,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])
    Str2 = paste(unique(strsplit(Str2,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])

    V(gg)[i]$TopOnto_OVG = as.character(Str1);
    V(gg)[i]$TopOnto_OVG_HDO_ID = as.character(Str2);


  }
  return(gg)
}
#Add topOnto_ov_P140papers
annotate_topOnto_ov_P140papers<-function(gg,par,dis){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"TopOnto_OV_PAPERS")
  gg <- removeVertexTerm(gg,"TopOnto_OV_PAPERS_HDO_ID")

  #--- Set Disease (geneRIF db) attributes in .gml graph
  set.vertex.attribute(gg,"TopOnto_OV_PAPERS",V(gg),"")
  set.vertex.attribute(gg,"TopOnto_OVG_PAPERS_HDO_ID",V(gg),"")

  #par    <- read.table("flatfile_human_gene2HDO.parentTerm.csv",sep="\t",skip=1,strip.white=T,quote="")
  #dis    <- read.table("flatfile_human_ov_PAPERS.csv",sep="\t",skip=1,strip.white=T,quote="")

  dis    <- rbind(dis,par)

  disIDS <- dis[,3]
  disn <-getDiseases()
  dtype<- getDType()

  for( i in 1:length(ids) ){

    ind1 = which(disIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){

      disv <- as.vector(dis[ind1,1]);

      indx <- match(disv,disn)

      for( j in 1:length(disv) ){

        if( !is.na(indx[j]) ){

          if( Str1 == "" ) { Str1 <- as.character(dtype[indx[j]]) }
          else {
            Str1 <- paste(c(Str1,as.character(dtype[indx[j]])),collapse=COLLAPSE[c]) }

          if( Str2 == "" ) { Str2 <- as.character(disn[indx[j]]) }
          else {
            Str2 <- paste(c(Str2,as.character(disn[indx[j]])),collapse=COLLAPSE[c]) }
        }

      }
    }

    Str1 = paste(unique(strsplit(Str1,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])
    Str2 = paste(unique(strsplit(Str2,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])

    V(gg)[i]$TopOnto_OV_PAPERS= as.character(Str1);
    V(gg)[i]$TopOnto_OV_PAPERS_HDO_ID = as.character(Str2);


  }
  return(gg)
}
#Add SCHanno synaptic functional groups
#' Title
#'
#' @param gg
#' @param anno
#'
#' @return
#' @export
#'
#' @examples
#' cid<-match('Presynaptic',getCompartments()$Name)
#' t<-getAllGenes4Compartment(cid)
#' gg<-buildFromSynaptomeByEntrez(t$HumanEntrez)
#' afile<-system.file("extdata", "SCH_flatfile.csv", package = "AnNet")
#' dis    <- read.table(afile,sep="\t",skip=1,header=F,strip.white=T,quote="")
#' #agg<-annotate_SCHanno(gg,dis)
annotate_SCHanno<-function(gg,anno){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"SCHanno")

  #--- Set Family attributes in .gml graph
  set.vertex.attribute(gg,"SCHanno",V(gg),"")

  #anno    <- read.table("SCH_flatfile.csv",sep="\t",skip=1,strip.white=T)
  annoIDS <- as.character(anno[,3])

  type <- unique(unlist(strsplit(as.character(unique(anno[,2])),",")))


  for( i in 1:length(ids) ){

    ind1 = which(annoIDS==ids[i])

    Str <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str <- as.character(anno[ind1[1],2]) }
      else { Str <- paste(as.character(anno[ind1,2]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$Schanno = as.character(Str);


  }
  return(gg)
}
#Add CHUA synaptic functional groups
annotate_CHUA<-function(gg,anno){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"chua")

  #--- Set Family attributes in .gml graph
  set.vertex.attribute(gg,"chua",V(gg),"")

  #anno    <- read.table("flatfile_chua.csv",sep="\t",skip=1,strip.white=T)
  annoIDS <- as.character(anno[,3])

  type <- unique(unlist(strsplit(as.character(unique(anno[,2])),",")))

  for( i in 1:length(ids) ){

    ind1 = which(annoIDS==ids[i])

    Str <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str <- as.character(anno[ind1[1],2]) }
      else { Str <- paste(as.character(anno[ind1,2]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$chua = as.character(Str);


  }
  return(gg)
}
#Add InterPro Family and Domain synaptic functional groups
annotate_Interpro<-function(gg,annoF,annoD){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"InterProFamilyID")
  gg <- removeVertexTerm(gg,"InterProFamily")

  gg <- removeVertexTerm(gg,"InterProDomainID")
  gg <- removeVertexTerm(gg,"InterProDomain")

  #--- Set interproFamily attributes in .gml graph
  set.vertex.attribute(gg,"InterProFamilyID",V(gg),"")
  set.vertex.attribute(gg,"InterProFamily",V(gg),"")
  set.vertex.attribute(gg,"InterProDomainID",V(gg),"")
  set.vertex.attribute(gg,"InterProDomain",V(gg),"")

  #annoF    <- read.table("flatfile.interpro.Family.csv",sep="\t",skip=1,strip.white=T)
  annoFIDS <- as.character(annoF[,3])

  typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))

  #annoD    <- read.table("flatfile.interpro.Domain.csv",sep="\t",skip=1,strip.white=T)
  annoDIDS <- as.character(annoD[,3])

  typeD <- unique(unlist(strsplit(as.character(unique(annoD[,2])),",")))


  for( i in 1:length(ids) ){

    ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$InterProFamilyID = as.character(Str2);
    V(gg)[i]$InterProFamily   = as.character(Str1);


    ind1 = which(annoDIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoD[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoD[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoD[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoD[ind1,1]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$InterProDomainID = as.character(Str2);
    V(gg)[i]$InterProDomain   = as.character(Str1);


  }

  return(gg)
}
#Add Core PSD and Pre-synpatic compartmental genes
annotate_compartments<-function(gg,preSet,psd95Set){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"COREPRE")
  gg <- removeVertexTerm(gg,"CORESPD")


  #---ADD VIP gene lists
  # VIP <- vector(length=2)
  # VIP[1] <- "SynsysBaits.csv"
  # VIP[2] <- "CorePSD95Complex.csv"
  # set1 <- read.table(sprintf("%s/%s",OUT[4],VIP[1]),sep="\t",header=F)[[1]]
  # set2 <- read.table(sprintf("%s/%s",OUT[4],VIP[2]),sep="\t",header=F)[[1]]
  set1 <- preSet
  set2 <- psd95Set

  set.vertex.attribute(gg,"COREPRE",V(gg),"")
  for( i in 1:length(ids) ){

    ind1 = which(set1==ids[i])

    Str <- "";

    if( length(ind1) != 0 ){
      Str <- "YES"
    }

    V(gg)[i]$COREPRE = as.character(Str);


  }

  set.vertex.attribute(gg,"CORESPD",V(gg),"")
  for( i in 1:length(ids) ){

    ind1 = which(set2==ids[i])

    Str <- "";

    if( length(ind1) != 0 ){
      Str <- "YES"
    }

    V(gg)[i]$COREPSD = as.character(Str);


  }
  return(gg)
}
#Add Bridgeness Regions for each algorithm
annotate_bridgeness_regions<-function(gg,str){
  ids = V(gg)$name
  #str <- sprintf("%s/%s/REGIONS/",OUT[4],subDIR[S])

  files <- list.files(str)

  fn <- gsub(".csv","",files)

  for( i in 1:length(files) ){

    gg <- removeVertexTerm(gg,fn[i])
    #gg <- removeVertexTerm(gg,gsub("_","",fn[i]))

    if( file.exists(sprintf("%s/%s",str,files[i])) ){

      ff <- read.table(sprintf("%s/%s",str,files[i]),sep="\t",header=F)

      gg <- set.vertex.attribute(gg,fn[i],V(gg),ff[match(ff[,1],ids),2])

    }
  }
  return(gg)
}
#Add GO MF
#' Title
#'
#' @param gg
#' @param annoF
#'
#' @return
#' @export
#'
#' @examples
annotate_go_mf<-function(gg,annoF){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"GO_MF")
  gg <- removeVertexTerm(gg,"GO_MF_ID")

  #--- Set Disease (geneRIF db) attributes in .gml graph
  set.vertex.attribute(gg,"GO_MF",V(gg),"")
  set.vertex.attribute(gg,"GO_MF_ID",V(gg),"")

#  annoF    <- read.table("flatfile.go.MF.csv",sep="\t",skip=1,strip.white=T,quote="")
  annoFIDS <- as.character(annoF[,3])

  typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))

  for( i in 1:length(ids) ){

    ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$GO_MF_ID = as.character(Str2);
    V(gg)[i]$GO_MF    = as.character(Str1);

  }
  return(gg)
}
#Add GO BP
#' Title
#'
#' @param gg
#' @param annoF
#'
#' @return
#' @export
#'
#' @examples
annotate_go_bp<-function(gg,annoF){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"GO_BP")
  gg <- removeVertexTerm(gg,"GO_BP_ID")

  #--- Set Disease (geneRIF db) attributes in .gml graph
  set.vertex.attribute(gg,"GO_BP",V(gg),"")
  set.vertex.attribute(gg,"GO_BP_ID",V(gg),"")

  #annoF    <- read.table("flatfile.go.BP.csv",sep="\t",skip=1,strip.white=T,quote="")
  annoFIDS <- as.character(annoF[,3])

  typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))

  for( i in 1:length(ids) ){

    ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$GO_BP_ID = as.character(Str2);
    V(gg)[i]$GO_BP    = as.character(Str1);

  }
  return(gg)
}
#Add GO CC
#' Title
#'
#' @param gg
#' @param annoF
#'
#' @return
#' @export
#'
#' @examples
annotate_go_cc<-function(gg,annoF){
  ids = V(gg)$name
  gg <- removeVertexTerm(gg,"GO_CC")
  gg <- removeVertexTerm(gg,"GO_CC_ID")

  #--- Set Disease (geneRIF db) attributes in .gml graph
  set.vertex.attribute(gg,"GO_CC",V(gg),"")
  set.vertex.attribute(gg,"GO_CC_ID",V(gg),"")

  #annoF    <- read.table("flatfile.go.CC.csv",sep="\t",skip=1,strip.white=T,quote="")
  annoFIDS <- as.character(annoF[,3])

  typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))

  for( i in 1:length(ids) ){

    ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";

    if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }

    }

    V(gg)[i]$GO_CC_ID = as.character(Str2);
    V(gg)[i]$GO_CC    = as.character(Str1);

  }
  return(gg)
}
#Add celltypes
annotate_celltypes<-function(gg,files){
  ids = V(gg)$name
  #files <- list.files("./")
  files <- files[grepl("celltypes_",files)]

  fn <- gsub(".csv","",files)
  fn <- gsub("celltypes_","",fn)
  fn <- sprintf("CellTypes_%s",fn)

  gg <- loopOverFiles(gg, files, fn, ids, FALSE)
  return(gg)
}
#Add pathways
annotate_pathways<-function(gg,files){
  ids = V(gg)$name
  #files <- list.files("./")
  files <- files[grepl("Pathways_",files)]

  fn <- gsub(".csv","",files)
  fn <- gsub("Pathways_","",fn)
  fn <- sprintf("PathWays_%s",fn)

  gg <- loopOverFiles(gg, files, fn, ids, TRUE)
  return(gg)
}
