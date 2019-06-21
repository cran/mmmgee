# 1 "independence",
# 2 "ar1", alpha.new
# 3 "exchangeable", alpha.new
# 4 "m-dependent", alpha.new, Mv=length(alpha.new)
# 5 "unstructured",  alpha.new (alpha.new ist vectorised upper tri part der Korr Matrix,
#			z.B. alpha.new= a,b,c, dann ist Korr=( 	1,a,b
#										a,1,c
# 										b,c,1)
# 6 "fixed", mat (the full korr matrix of maximal size)
# 7 "userdefined", alpha.new, struct.vec (which is the vectorised representation of the upper triangle of the user provided corr.mat structure


getAlphaInvGeneral<-function(alpha.new, row.vec, col.vec, len, includedvec, inclsplit, cor.match, struct.vec=NULL) {
	K<-length(len)
	###
	#this part should in the end be at the start of the main function, it needs to be computed
	#only once
	matrixID<-vector(length=K,mode="numeric")
	matrixIDset<-NULL
	matrixIDlen<-NULL
	matrixIDincl<-list()
	zz<-0
	for(i in 1:K) {
		matrixID[i]<-as.numeric(paste(c(len[i],as.numeric(inclsplit[[i]])),collapse=""))
		if(!(matrixID[i]%in%matrixIDset)) {
			zz<-zz+1
			matrixIDset<-c(matrixIDset,matrixID[i])
			matrixIDlen<-c(matrixIDlen,len[i])
			matrixIDincl[[zz]]<-inclsplit[[i]]
		} 
	}
	NMat<-length(matrixIDset)
	#biggest possible correlation matrix
	unstr.row <- NULL
	unstr.col <- NULL
	ml <- max(len)
	for (i in 2:ml) {
	      unstr.row <- c(unstr.row, 1:(i - 1))
      	unstr.col <- c(unstr.col, rep(i, each = i - 1))
	}
	unstr.row <- c(unstr.row, 1:ml)
	unstr.col <- c(unstr.col, 1:ml)
	#
	###


	if(cor.match == 2) { #AR
		xvec<-alpha.new^(unstr.col-unstr.row)
	}

	if(cor.match == 3) { #exchangeable
		xvec<-c( rep(alpha.new,ml*(ml-1)/2) , rep(1, ml))
	}

	if(cor.match == 4) { #m-dependent, with m=Mv
		xveclen<-ml*(ml-1)/2
		xvec<-rep(0,xveclen)
		delta<-(unstr.col-unstr.row)[1:xveclen] #the remaining entries are the diagonal indices, must be excluded here
		for(j in 1:length(alpha.new)) xvec[delta==j]<-alpha.new[j]
		xvec<-c(xvec, rep(1, ml))
	}

	if(cor.match == 5) { #unstructured
		xvec <- c(alpha.new, rep(1, ml))
	}

	if(cor.match == 6) { #fixed
		#in this case, alpha.new must be prepared in the main function as vectorized upper triangle matrix of 
		#corr.mat
		xvec <- c(alpha.new, rep(1, ml))
	}
	
	if(cor.match == 7) { #user defined
		xvec <- rep.int(0, length(struct.vec))
		for (i in 1:length(alpha.new)) {
			xvec[which(struct.vec == i)] <- alpha.new[i]
		}
		xvec <- c(xvec, rep(1, ml))
	}

	biggestMat <- forceSymmetric(sparseMatrix(i = unstr.row, j = unstr.col, x = xvec))
	#Inverse berechnen
	mat.inverses <- vector(mode="list",length=NMat)
	for(i in 1:NMat) {
		zeros<-!matrixIDincl[[i]]
		anyzeros<-any(zeros)
		#in case of EXCH and AR1, both without zeros, more efficient algorithms then solve may be used:
		#they are implemented for matrix dimension > 2, for dim 1 and 2 solve is sufficient anyways.
		if(cor.match == 3 & !anyzeros & matrixIDlen[i]>2) {
			mat.inverses[[i]]<-solveEXCH.Vek(alpha.new,matrixIDlen[i])
		} else if(cor.match == 2 & !anyzeros & matrixIDlen[i]>2) {
			mat.inverses[[i]]<-solveAR1.Vek(alpha.new,matrixIDlen[i])
		} else {
			M <- biggestMat[1:matrixIDlen[i], 1:matrixIDlen[i]]
			#The following step is required when the data just have missings and we do NOT account for that using residual weights.
			#If, however, residual weights are used, this means for the remaining data that implicitly there are other obervations,
			#and thus the weighting matrix resulting from the working correlation should account for these observations.
			#See Preisser, J. S., Lohman, K. K., & Rathouz, P. J. (2002). Performance of weighted estimating equations for longitudinal binary data with drop-outs missing at random. Statistics in medicine, 21(20), 3035-3054.
			if(anyzeros) { #& !residualWeightsInUse
				M[zeros,]<-0 
				M[,zeros]<-0
				diag(M)<-1
				#This is very important. The correct inverse needs to be calcualted from a correlation matrix covering only the non-missing values.
				#This is not the same as to invert and then remove the missing values, because entries in the inverse are potentially connected to
				#all entries in the original matrix. Therefore: First remove missing values, then calculate inverse. With dummy rows, this is technically the same
				#as setting all correlations that involve missing values to zero, as is done in the code here, and then inverting.
				#das ist sehr wichtig, denn die Inverse enthaelt sonst Eintraege, die aufgrund eigentlich fehlender Werte verwendet werden
				#und dadurch haben auch die nicht fehlenden ein anderes Gewicht als sie haetten, wenn man nur mit den vorhandenen die 
				#Korrelationsmatrix bildet und invertiert
			}
			mat.inverses[[i]]<-as.vector(solve(M))
		}
		#if(i%in%c(1,2,3)) print(mat.inverses[[i]])
	}

	mat.finder <- match(matrixID, matrixIDset)
	corr.vec <- unlist(mat.inverses[mat.finder])
 	#Rinverse<-as(sparseMatrix(i = row.vec, j = col.vec, x = corr.vec), "symmetricMatrix") 
		#this may not work, if due to numeric inaccuracies the corr.vec values are not exactly the same at symmetric positions.
		#which may be the case for large blocks, though the inaccuracies are in the order of the machine accuracy
		#better use forceSymmetric
	Rinverse<-forceSymmetric(sparseMatrix(i = row.vec, j = col.vec, x = corr.vec))
	return(list(Rinverse,biggestMat))
}


#for mdim>2
solveEXCH.Vek<-function(rho,mdim) {
	offdiagonal<- -rho/((1-rho)*(1+(mdim-1)*rho))
	diagonal<- (1+(mdim-2)*rho)/((1-rho)*(1+(mdim-1)*rho))
	return(c( rep(c(diagonal,rep(offdiagonal,mdim)),mdim-1),diagonal))
}

#for mdim>2
solveAR1.Vek<-function(rho,mdim) {
	a<-1/(1-rho^2)
	b<- -rho*a
	c<- (1+rho^2)*a
	zerofill<-rep(0,mdim-2)
	return(c(a, rep(c(b,zerofill,b,c),mdim-2) ,c(b,zerofill,b,a) ))
}




