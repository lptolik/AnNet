library(igraph);
library(lattice);
library(methods);
library(poweRlaw);
library(scales);
library(grid);
library(latex2exp);
library(stringr);

##---WIDTH and HEIGHT for plots
WIDTH=480
HEIGHT=480

##set igraph as S4
setClass("poweRlaw")

## Gene Set Analysis (GSA) object
setClass(Class="law",representation(
  fit="displ",
  p="numeric",
  alpha="numeric",
  SDxmin="numeric",
  SDalpha="numeric")
)

##log sequence of numbers
lseqBy <- function(from=1, to=100000, by=1, length.out=log10(to/from)+1) {
  tmp <- exp(seq(log(from), log(to), length.out = length.out))
  tmp[seq(1, length(tmp), by)]
}

changeSciNot <- function(n) {

  n <- format(n, scientific = TRUE)

  oo <- strsplit(as.character(n),"e")

  out <- vector(length=length(oo))

  out[1] <- TeX(sprintf("10^{%s}",sub("-0?","-",oo[[1]][2])))

  for( i in 2:length(oo) ){

    if(grepl("-",oo[[i]][2])){
      out[i] <- TeX(sprintf("10^{%s}",sub("-0?","-",oo[[i]][2])))
    } else {
      out[i] <- TeX(sprintf("10^{%s}",sub("\\+0?","",oo[[i]][2])))
    }

  }

  return(out)


}

#' Fit Power Law to degree distribution.
#'
#' @param DEG degree distribution
#' @param title
#' @param Nsim
#' @param DATAleg
#' @param WIDTH
#' @param HEIGHT
#'
#' @return
#' @export
#'
#' @examples
#' ##No: of bootstrap iterations
#' bootStrap    <- vector(length=3)
#' bootStrap[1] <- 100
#' bootStrap[2] <- 1000
#' bootStrap[3] <- 5000
#'
#' b=2 #we'll stick with 1000 iteteractions
#'
#' ##Legend Titles
#' Legend <- vector(length=2)
#' Legend[1] <- "Presynaptic PPI"
#' Legend[2] <- "PSP PPI"
#'
#' l=2;
#'
#' gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[3],subDIR[S],subDIR[S]),format="gml")
#' pFit <- FitDegree( as.vector(igraph::degree(graph=gg)), subDIR[S], bootStrap[b], Legend[l], WIDTH, HEIGHT )

FitDegree <- function(DEG, title, Nsim, DATAleg, WIDTH, HEIGHT ){

  #WIDTH=480
  #HEIGHT=480

  #tmp = max(DEG)
  #Max = 4*tmp

  DEG <- DEG[DEG > 0]

  data <- DEG

  m_pl = displ$new(data)

  est = estimate_xmin(m_pl)

  m_pl$setXmin(est)

  gof <- bootstrap_p(m_pl, no_of_sims = Nsim, threads=4)

  x_lab="k"    ##degree
  y_lab="P(k)" ## the CDFs P(k) for the PPI network data
  leg_x = max(data)
  leg_y = 1

  png(filename=sprintf("PLOTS/%s_cdf.png",title), width=WIDTH, height=HEIGHT, units="px")

  d = plot(m_pl,draw=F)

  Xmax <- max(d$x) - max(d$x)*0.5

  ##build y-axis labels
  yTICKS  = round(lseqBy(min(d$y),1,0.5),4)
  yLABELS = changeSciNot(yTICKS)

  plot(m_pl, xlab=sprintf("%s",x_lab), ylab=y_lab,
       panel.first=grid(col="grey60"),
       pch=22, bg='black', axes = F, cex.lab = 1.5, yaxt='n' )
  box(col='black')

  axis(1, cex.axis = 1.5, font = 1.5, family = 'arial')
  axis(2, cex.axis = 1.5, font = 1.5, family = 'arial', at=yTICKS, labels=yLABELS)

  lines(m_pl, col=2, lwd=3)

  S1 = round(m_pl$xmin,2)
  S2 = round(m_pl$pars,2)
  S3 = round(gof$p,2)

  sdS1 = round(sd(gof$bootstraps$xmin),0)
  sdS2 = round(sd(gof$bootstraps$pars),2)

  errS1 = str_sub(as.character(sdS1),-1,-1)
  errS2 = str_sub(as.character(sdS2),-1,-1)


  fitl <- TeX(sprintf("Power-law $\\alpha = %.2f(%s), $k_{min} = %.0f(%s)",S2,errS2,S1,errS1))

  legend("bottomleft",c(DATAleg,fitl),lty=c(1,1),lwd=c(4,4),col=c('black',2),merge=TRUE, cex = 1.5)

  #legend("topright",c(DATAleg,expression(paste(P(k), " = ",frac(k^alpha,sigma1(alpha,k[min]))))),lty=c(1,1),lwd=c(4,4),col=c('black',2),merge=TRUE, cex = 1.5)
  #text(x=Xmax,y=0.2,substitute(paste(k[min], " = ", s1, ", ", alpha, " = ", s2), list(s1=S1,s2=S2)),cex=0.9)

  dev.off()


  return(new("law",
             fit=m_pl,
             p=as.numeric(gof$p),
             alpha=as.numeric(est$pars),
             SDxmin=as.numeric(sd(gof$bootstraps$xmin)),
             SDalpha=as.numeric(sd(gof$bootstraps$pars))
  )
  )


}
