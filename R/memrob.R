
#' function to get the member robustness from the consensus matrices #BUG fixed 14/03/12 TIS
#'
#' @param x
#' @param rm
#'
#' @return
#' @import clusterCons
#'
#' @examples
memrob <- function(x,rm=data.frame()){
  if(class(x) == 'consmatrix'){
		cmref <- x@rm
	}
	else{
		if(length(rm)==0){stop('You need to specify a reference matrix for a merge consensus matrix')}
		else{
			cmref <- rm
		}
	}

	consensus <- x@cm

  #BUG - fixed to ensure deals with >100 clusters
	mem_rob = matrix(0,dim(consensus)[1],length(levels(as.factor(cmref$cm))),dimnames = list(row.names(consensus),1:length(levels(as.factor(cmref$cm)))))

	for(k in 1:length(levels(as.factor(cmref$cm)))){ #BUG - fixed to ensure deals with >100 clusters
		for(i in 1:dim(consensus)[1]){
			Ik = row.names(cmref)[cmref$cm==k] #where k is the cluster number

			ind = Ik[Ik != row.names(consensus)[i]] # exclude the index for i = j if it is there

			sigma = apply(as.matrix(consensus[ind,i]),2,sum) #perform the sigma sum on index

			ei = row.names(consensus)[i] # get the current member we are checking

			Nk = summary(as.factor(cmref$cm),maxsum=10000)[k] # get the current cluster size note this is limited to a max of 10000 custers


			if(sum(ei == Ik) ==1){		#if ei is a member of Ik
				mik = (1/(Nk-1))*sigma
			}
			else{				#if ei is not a member of Ik
				mik = (1/Nk)*sigma
			}
			mem_rob[i,k] = mik
		}
	}
	#what you might want to do here is have the full object in a slot and output the mem_rob for the ref clustering (which is what you actually want)
	mem_rob_list <- list();
	for(i in 1:dim(mem_rob)[2]){
		cl_mem_rob <- (mem_rob[(cmref==i),i])
		current_list <- data.frame(sort(cl_mem_rob,dec=TRUE))
		names(current_list) <- 'mem_rob'
		current_mem_rob_list <- new('memroblist',mrl=current_list);
		mem_rob_list[paste('cluster',i,sep='')] <- current_mem_rob_list
	}
	mem_rob_list['resultmatrix']<- new('memrobmatrix',mrm=mem_rob);
	mem_rob_list['algo']<- x@a;
	mem_rob_list['type']<- class(x);
	return(mem_rob_list)
}

#' Calculate cluster robustness from consensus matrix
#'
#' @param gg igroph object
#' @param alg clustering algorithm
#' @param conmat consensus matrix
#'
#' @return
#' @export
#'
#' @examples
getRobustness<-function(gg,alg,conmat){

  if(!alg%in%igraph::vertex_attr_names(gg)){
    stop(paste("Membership for the ",alg,
               "algorithm should be stored in the graph.\n",
               "See calcClustering.\n"))
  }
  rm<-data.frame(cm=as.numeric(igraph::get.vertex.attribute(gg,alg,V(gg))))
  cm           <- data.frame(conmat);
  names(cm)    <- rownames(rm);
  rownames(cm) <- rownames(rm);
  cm           <- as.matrix(cm);
  kk     = max(as.numeric(as.vector(rm$cm)));
  ##--- make the consensus matrix object for clusterCons so you can use its functions
  out <- new('consmatrix', cm=cm, rm=rm, k=kk, a=alg);

  ##--- get cluster robustness values from clusterCons
  cr <- clusterCons::clrob(out);

  ##--- the scaled cluster robustness values
  crScales <- cr$rob
  crScales <- (crScales-min(crScales))/(max(crScales)-min(crScales))
  oo<-data.frame(alg=as.character(rep(alg,length(rownames(cr)))),
             C=as.numeric(rownames(cr)),
             Cn=as.numeric(table(rm[,1])),
             Crob=as.numeric(cr$rob),
             CrobScaled=as.numeric(crScales))

}
