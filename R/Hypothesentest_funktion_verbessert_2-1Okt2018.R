#' @title Hypothesis Tests for Linear Contrasts in Multiple Marginal GEE Models
#' @description Global hypothesis tests, multiple testing procedures and simultaneous confidence intervals for multiple linear contrasts of regression
#'    coefficients in a single generalized estimating equation (GEE) model or across multiple GEE models.
#' @param x a \code{geem2} object fitted with \code{geem2} or a list of \code{geem2}. In the latter case, the \code{geem2} objects must be different
#'	models calculated with data from the same subjects. In particular, the parameter \code{id} in the call to \code{geem2} must refer to the same subjects in each model. 
#' @param L a contrast matrix defining a contrast for the stacked vector of regression coefficients of the marginal models, or a list of contrast matrices.
#'	In the latter case, the list must contain one matrix for each model listed in \code{x}, in the same order as the models. When using the the score test \code{x} must it be a list.
#' @param r right hand side vector of the null hypothesis or a list of vectors resembling the right hand side of the null hypothesis. If not specified \code{r=0} is assumed. See details.
#' @param statistic either \code{"wald"} or \code{"score"}, see details. The default is \code{"wald"}.
#' @param type either \code{"maximum"} or \code{"quadratic"}, see details. The default is \code{"maximum"}.
#' @param asymptotic logical, if \code{TRUE} the reference distribution for the maximum-type Wald test statistic is a multivariate normal distribution and the reference distribution for 
#' 	the quadratic form Wald test statistic is a chi-squared distribution. If \code{FALSE}, a multivariate t-distribution or an F-distribution is used instead.
#'	Ignored for the Score test, see details.
#' @param biascorr logical indicating whether the Mancl and DeRouen Bias correction should be used when estimating the joint covariance matrix via \code{\link{mmmgee}}.
#' @param closed.test logical, if \code{TRUE}, multiplicity adjusted p-values based on a closed test procedure using the selected type of test are calculated. With \eqn{k} hypotheses this 
#'	involves the computation of \eqn{2^k} tests, which may require considerable computation time.
#' @param conf.int logical. If \code{TRUE} simultaneous confidence intervals corresponding to a single step maximum-type test are calculated using a multivariate
#'	normal or t approximation, depending on \code{asymptotic}.
#' @param conf.level the nominal simultaneous coverage probability of the confidence intervals.
#' @param alternative one of \code{"undirected"}, \code{"greater"}, or \code{"less"}. Determines the direction of maximum-type tests and of confidence intervals.
#'	The default is \code{"undirected"}.
#' @param denomDF Defaults to \code{NULL}. In that case, denominator degrees of freedom for the multavariate t-distribution or F-distribution are calculated as \code{min(n-p)},
#'	where \code{n} and \code{p} are vectors of the number of independent clusters and the number of regression coefficients in the models in \code{x}.
#'	Alternatively, a numeric value may be entered to be used as denominator degrees of freedom. 
#' @param scaled.F logical. If \code{TRUE} and \code{type="quadratic"} and \code{asymptotic=FALSE} a scaled F distribution similar as for Hotelling's test is used. Ignored otherwise.
#' @param tol tolerance limit for the convergence criterion to be passed to \code{\link{geem2}}. Only required when using the score test, where the models are refitted under the
#' 	restriction of the null hyptothesis.
#' @param ... additional arguments that are passed to \code{\link[mvtnorm]{pmvnorm}}, \code{\link[mvtnorm]{qmvnorm}}, \code{\link[mvtnorm]{pmvt}} and \code{\link[mvtnorm]{qmvt}}.
#'	In particular the algorithm to solve the multivariate normal or t-distribution integrals may be selected.
#'
#' @details The null hypothesis is \eqn{H0:L\beta=r} where \code{L} is a contrast matrix, \eqn{\beta}
#'	the stacked vector of regression coefficients rom the marginal models and \code{r} a real values right hand side vector.
#'	\code{L} can be specified as matrix or, if it is a block diagonal matrix with each block corresponding to 
#'	a contrast for one marginal GEE model, as list of the matrices on the diagonal. The right hand side \code{r} can be speficied as vector or as list of 
#'	vectors each corresponding to the part of the right hand side vector for one model.
#'
#' 	When choosing \code{statistic="wald"} and \code{type="maximum"}, the maximum of the standardized entries of \eqn{L\hat{\beta}} is used as test statistic and the
#'	p-value is calculated from a multivariate normal or t-distribution (depending on \code{asymptotic} being \code{TRUE} or \code{FALSE}) with correlation matrix
#'	estimated for \eqn{L\hat{\beta}}. For the t-distribution, denominator degrees of freedom are used as specified in \code{denomDF}.
#' 	When choosing \code{statistic="wald"} and \code{type="quadratic"}, a quadratic form of \eqn{L\hat{\beta}} and the inverse of the estimated covariance matrix
#'	of \eqn{L\hat{\beta}} is used as test statistic and the p-value is calculated from a 
#'	chi-squared distribution or an F-distribution (depending on \code{asymptotic} being \code{TRUE} or \code{FALSE}).
#'
#'	With \code{statistic="score"}, generalized score tests are calculated by replacing \eqn{L\hat{\beta}} by the first order approximation \eqn{LAU} where \eqn{U}
#'	is the stacked estimating equation (the score) and \eqn{A} is the negative inverse of the matrix
#'	of first derivatives of \eqn{U}, both evaluated at the location of constrained estimates for \eqn{\beta} under the null hypothesis. 
#'	Analogous to the Wald statistic, a maximum-type and a quadratic form score test are available. For the score test the option \code{asymptotic} is ignored
#'	and the reference distribution is multivariate normal or chi-squared.
#'
#' @note Calculating the generalized score test requires refitting the models under the constraint of the null hypothesis.
#'	The function \code{\link[stats]{update}} is used for this task. It will use the function calls as stated in the slot \code{call} of the \code{geem2} objects.
#'	There is one important point to notice: The function \code{update} will first look for any components of the fitted object in the environment from which it was
#'	called, which is an internal function of the package.
#'	Within this internal function the variables 'Modelle', 'L.list', 'r.list', 'tol', 'type', 'alternative' and 'biascorr' are used. If any component of the model
#'	object happens to have one of these names, (e.g. if your data frame is called 'tol'), \code{update} will erronously try to use the internal object of that name.
#'
#'	A single value for the denominator degrees of freedom is calculated for the covariance matrix estiamate across all contrasts. In the close testing procedure,
#'   	this value is used for the degrees of freedom
#'	associated with the covariance matrix of any subset of contrasts.
#'
#'	Usual linear models or generalized linear models can be regarded as special case of GEE models and can be included in the analysis framework.
#'
#' @return A list with class \code{mmmgeetest} containing the following components, if required:
#' \describe{
#'	\item{\code{test}}{Contains a data frame with the test statistic, degrees of freedom (depending in the type of test) and the p-value. If closed test was required,
#'		a further data frame is reported with estimates, right hand side vector, unadjusted p-values and adjusted p-values for each line of \eqn{H0:L\beta-r=0}.}
#'	\item{\code{hypothesis}}{A list containing the contrast matrix \eqn{L} and the right hand side vector \eqn{r}.}
#'	\item{\code{conf.int}}{The simultaneous confidence intervals.}
#'	\item{\code{denomDF}}{The type and value of the denominator degrees of freedom used in the procedure.}
#'	\item{\code{mmmgee}}{The \code{mmmgee} object containing in particular the estimated covariance matrix for the coefficents of the models in \code{x}. See \code{\link{mmmgee}}.}
#' }
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @references Dennis D. Boos. On generalized score tests. The American Statistician, 1992, 46(4):327-333.
#' @references Lloyd A. Mancl, Timothy A. DeRouen. A covariance estimator for GEE with improved small sample properties. Biometrics, 2001, 57(1):126-134.
#'
#' @seealso \code{\link{geem2}}, \code{\link{mmmgee}}
#'
#' @examples
#' data(keratosis)
#' m1<-geem2(clearance~trt,id=id,data=keratosis,family=binomial,corstr="independence")
#' m2<-geem2(pain~trt,id=id,data=keratosis[keratosis$lesion==1,],family=gaussian,corstr="independence")
#' L1<-L2<-diag(1,4)[-1,]
#' mmmgee.test(x=m1,L=list(L1),statistic="wald",type="maximum")
#' mmmgee.test(x=m1,L=list(L1),statistic="score",type="maximum")
#' mmmgee.test(x=list(m1,m2),L=list(L1,L2),type="maximum",asymptotic=FALSE,biascorr=TRUE)
#' mmmgee.test(x=list(m1,m2),L=list(L1,L2),type="maximum",closed.test=TRUE)
#' mmmgee.test(x=list(m1,m2),L=list(L1,L2),type="maximum",asymptotic=FALSE,
#' 	alternative="less",conf.int=TRUE,denomDF=40)
#' mmmgee.test(x=list(m1,m2),L=list(L1,L2),type="quadratic",asymptotic=TRUE)
#' mmmgee.test(x=list(m1,m2),L=list(L1,L2),statistic="score",type="quadratic")
#' mmmgee.test(x=list(m1,m2),L=list(L1,L2),statistic="score",type="maximum")
#'
#' @export
mmmgee.test<-function(x,L=NULL,r=NULL,statistic=c("wald","score"),type=c("maximum","quadratic"),asymptotic=TRUE,biascorr=FALSE,closed.test=FALSE,conf.int=FALSE,conf.level=0.95,alternative=c("undirected","greater","less"),denomDF=NULL,scaled.F=FALSE,tol=10^(-8),...) {
	if(is.null(L)) {
		L<-vector(mode="list",length=length(x))
		for(i in 1:length(x)) {
			L[[i]]<-diag(1,length(x[[i]]$beta))
			rownames(L[[i]])<-paste("model",i,x[[i]]$coefnames,sep="_")
		}
	}
	if(class(x)=="geem2") {
		x<-list(x)
	} else {
		if(!is.list(x)) stop("x must be a geem2 object or a list of geem2 objects.")
		for(i in 1:length(x)) if(class(x[[i]])!="geem2") stop("x must be a list of geem2 objects fitted with geem2.")
	}

	statistic<-match.arg(tolower(statistic),choices=c("wald","score"))
	if(statistic=="wald") statistic<-"Wald"
	if(statistic=="score") statistic<-"Score"
	type<-match.arg(tolower(type),choices=c("maximum","quadratic"))
	#Set reference distribution:
	#approximation<-match.arg(approximation,choices=c("normal","t","Chisq","F","scaled.F")) ### Chi2,F, scaled.F auch dazu??
	#approximation=c("normal","t","Chisq","F","scaled.F")
	if(type=="maximum") {
		if(statistic=="Wald" & asymptotic==TRUE) approximation<-"normal"
		if(statistic=="Wald" & asymptotic==FALSE) approximation<-"t"
		if(statistic=="Score") approximation<-"normal"
	}
	if(type=="quadratic") {
		if(statistic=="Wald" & asymptotic==TRUE) approximation<-"Chisq"
		if(statistic=="Wald" & asymptotic==FALSE) approximation<-"F"
		if(statistic=="Wald" & asymptotic==FALSE & scaled.F==TRUE) approximation<-"scaled.F"
		if(statistic=="Score") approximation<-"Chisq"
	}
	#conf.int<-match.arg(conf.int,c("none","MaxZ","MaxT"))
	alternative<-match.arg(alternative,c("undirected","greater","less"))
	if(type=="quadratic" & alternative!="undirected") warning('alternative different from "undirected" is ignored for quadratic form tests.')
	if(statistic=="Score" & !is.list(L)) stop("L must be a list if using the Score test.")
	if(statistic=="Score" & !is.null(r)) if(!is.list(r)) stop("r must be a list if using the Score test.")

	if(is.list(L) & is.list(r) & length(L)!=length(r) ) stop("If lists, L and r must be of same length.")
	if(is.list(L)) LL<-bdiag(L) else LL<-L
	df1<-dim(LL)[1]
	if(is.list(r)) rr<-unlist(r) else rr<-r
	if(is.null(r)) rr<-rep(0,df1) 
	if(dim(LL)[1]!=length(rr)) stop("Dimensions of contrast matrix and right-hand side vector do not match.")
	
	if(is.list(L)) contrastnames<-unlist(lapply(L,rownames))
	if(is.matrix(L)) contrastnames<-rownames(L)
	if(is.null(contrastnames)) contrastnames<-1:df1

	#check ranks and if L%*%beta=r has a solution
	rankL<-qr(as.matrix(LL))$rank
	rankLr<-qr(as.matrix(cbind(LL,rr)))$rank
	if(rankL!=rankLr) stop("The system Lbeta = r has no solution.")
	if(type=="quadratic" & rankL!=dim(LL)[1]) stop("For quadratic form tests the contrast matrix must have full rank.")

	#if(!(denomDF=="minnp" | (is.numeric(denomDF) & denomDF>0) | denomDF=="satt")) stop('denomDF must be one of "minnp", "satt" or a positive number.') 
	xx<-mmmgee(x,biascorr=biascorr)
	if(dim(LL)[2]!=length(xx$beta)) stop("Dimensions of contrast matrix and vector of coeffiecients do not match.")
	
	est<-as.numeric(LL%*%xx$beta) #always included in output

	
	#set DF
	#Fuer mvt muss df2 ein Integer sein, in den Funktionen, die mvt verwenden, wird df2 abgerundet.
	if(is.null(denomDF)) {
		df2<-min(xx$n-xx$p)
		denomDFtype<-"minnp"
	} else {
		df2<-denomDF
		denomDFtype<-"user"
		if(statistic!="Wald" | approximation=="normal") warning("denomDF is only considered for the Wald test statistic with t or F approximation.")	
	}
	#if(denomDF=="minnp") {
	#	df2<-min(xx$n-xx$p)
	#	denomDFtype<-"minnp"
	#} else {
	#	if(denomDF=="satt") {
	#		df2<-xx$satt.df
	#		denomDFtype<-"satt"
	#	}
	#	if(is.numeric(denomDF)) {
	#		df2<-denomDF
	#		denomDFtype<-"user"
	#	}
	#}

	
		
	if(statistic=="Score") {
		global<-mmmscore.test(Modelle=x,L.list=L,r.list=r,tol=tol,type=type,alternative=alternative,biascorr=biascorr)
		Ldims<-sapply(L,function(x) dim(x)[1]) #wird fuer closed test benoetigt
	}
	if(statistic=="Wald") {
		if(type=="maximum") global<-mmmmax.test.simple(x=xx,L=LL,r=rr,alternative=alternative,approximation=approximation,df2=df2)
		if(type=="quadratic") global<-mmmchisq.test.simple(x=xx,L=LL,r=rr,approximation=approximation,df2=df2)
	
		#if(type=="Chisq" | type=="F" | type=="scaled.F") global<-mmmchisq.test.simple(x=xx,L=LL,r=rr,type=type,df2=df2)
		#if(type=="MaxT") global<-mmmmax.test.simple(x=xx,L=LL,r=rr,df2=df2,alternative=alternative)
		#if(type=="MaxZ") global<-mmmmax.test.simple(x=xx,L=LL,r=rr,df2=NULL,alternative=alternative)
	}

	#confint
	if(conf.int==TRUE) {
		alpha<-1-conf.level
		tail<-c("both.tails","lower.tail","upper.tail")[alternative==c("undirected","greater","less")]
		Vcontr<-LL%*%xx$V%*%t(LL)
		if(approximation%in%c("normal","Chisq")) {
			Q<-qmvnorm(p=1-alpha,tail = tail, sigma = cov2cor(as.matrix(Vcontr)),...)$quantile
			KIapprox<-"normal"
		}
		if(approximation%in%c("t","F","scaled.F")) {
			Q<-qmvt(p=1-alpha,tail = tail, sigma = cov2cor(as.matrix(Vcontr),...),df=floor(df2))$quantile
			KIapprox<-"t"
		}
		SE<-sqrt(diag(Vcontr))
		
		delta<-SE*Q
		if(alternative=="undirected") KI<-data.frame(contrast=contrastnames,estimate=est,low=est-delta,up=est+delta)
		if(alternative=="less") KI<-data.frame(contrast=contrastnames,estimate=est,low=-Inf,up=est+delta)
		if(alternative=="greater") KI<-data.frame(contrast=contrastnames,estimate=est,low=est-delta,up=Inf)
		#KI
		#Q
		#mc<-glht(model=xx,linfct=as.matrix(LL),df=df2)
		#confint(mc)
	} else {
		KI<-NULL
	}

	#closed test
	if(closed.test==TRUE) {
		S<-namen.fn(df1)[-1,,drop=FALSE] == 1
		S<-S[order(rowSums(S)),,drop=FALSE]
		#S[rowSums(S)==1,]
		k<-dim(S)[1]
		inClosedTest<-rep(TRUE,k)
		inClosedTest[k]<-TRUE #global test is always included
		p.intersect<-rep(NA,k)
		
		p.intersect[k]<-global$test$p

		#layer<-1+rowSums(!S) #there are as many layers as dimensions (contrasts)
		#with k we get the full table and hence the global intersection test
		#i<-20
		if(k>1) {
		    for(i in (k-1):1) {
			#Remove redundant hypotheses like in a Shaffer type test:
			#Check if the hypothesis has been tested before
			L.temp<-LL[S[i,],,drop=FALSE]
			stopHcheck<- rankL==dim(LL)[1] #if rank is full, we do not need to check for redundant hypotheses
			testdone<-FALSE
			i.prev<-k
			while(!stopHcheck) {
			#for(i.prev in k:(i+1)) {
				L.prev<-LL[S[i.prev,],,drop=FALSE]
				if(qr(as.matrix(L.temp))$rank==qr(as.matrix(rbind(L.prev,L.temp)))$rank) {
					#then the hypothesis is contained in a previous hypothesis
					p.intersect[i]<-p.intersect[i.prev]
					testdone<-TRUE
					inClosedTest[i]<-FALSE
					stopHcheck<-TRUE
				}
				i.prev<-i.prev-1
				if(i.prev==i) stopHcheck<-TRUE
			}	

			#
			if(testdone==FALSE) {
				select<-S[i,]

				if(statistic=="Score") {
					L.temp<-L
					x.temp<-x
					r.temp<-r
					for(j in length(L):1) { #backwards because we remove models from x by indexing and the index would change if we start in the front
						start<-c(0,cumsum(Ldims))[j]+1
						stop<-cumsum(Ldims)[j]
						subselect<-select[start:stop]
						L.temp[[j]]<-L.temp[[j]][subselect,,drop=FALSE]
						r.temp[[j]]<-r.temp[[j]][subselect]
						if(dim(L.temp[[j]])[1]==0) {
							x.temp[[j]]<-NULL
							L.temp[[j]]<-NULL
						}
					}
					#L.temp
					test.intersect<-mmmscore.test(Modelle=x.temp,L.list=L.temp,r.list=r.temp,tol=tol,type=type,alternative=alternative,biascorr=biascorr)
				}

				if(statistic=="Wald") {
					#L.temp<-LL[S[i,],,drop=FALSE] #steht jetzt oben
					r.temp<-rr[S[i,]]
					#if(type=="Chisq" | type=="F" | type=="scaled.F") test.intersect<-mmmchisq.test.simple(x=xx,L=L.temp,r=r.temp,type=type,df2=df2)
					#if(type=="MaxT") test.intersect<-mmmmax.test.simple(x=xx,L=L.temp,r=r.temp,df2=df2,alternative=alternative)
					#if(type=="MaxZ") test.intersect<-mmmmax.test.simple(x=xx,L=L.temp,r=r.temp,df2=NULL,alternative=alternative)
					if(type=="maximum") test.intersect<-mmmmax.test.simple(x=xx,L=L.temp,r=r.temp,alternative=alternative,approximation=approximation,df2=df2)
					if(type=="quadratic") test.intersect<-mmmchisq.test.simple(x=xx,L=L.temp,r=r.temp,approximation=approximation,df2=df2)
				}

				#test.intersect<-test.fun(x)
				p.intersect[i]<-test.intersect$test$p
			    }
			}
		}
		#unadjusted p-values
		p<-p.intersect[rowSums(S)==1]
		#adjusted p-values
		p.adj<-rep(NA,df1)
		for(i in 1:df1) p.adj[i]<-max(p.intersect[S[,i]])
		closed.test.out=data.frame(contrast=contrastnames,estimate=est,rhs=rr,p.unadj=p,p.adjusted=p.adj)
		redundancies<-data.frame(S,redundant=!inClosedTest)
		names(redundancies)<-c(paste("Contrast",1:df1),"redundant")
	} else {
		closed.test.out<-NULL
		redundancies<-NULL
	}

	#Output
	test.out=list(
		global=global$test,
		closed.test=closed.test.out,
		statistic=statistic,type=type,approximation=approximation,
		alternative=ifelse(type=="maximum",alternative,"undirected") )
		
	if(conf.int==TRUE) {
		conf.int.out=list(conf.int=KI,conf.level=conf.level,quantile=Q,alternative=alternative,approximation=KIapprox)
	} else {
		conf.int.out<-NULL
	}
	denomDF.out<-list(df2=df2,denomDFtype=denomDFtype)
	LL.out<-as.matrix(LL)
	rownames(LL.out)<-contrastnames
	#out<-list(test=test.out,hypothesis=list(contr.matrix=LL.out,rhs=rr,redundancies=redundancies),conf.int=conf.int.out,denomDF=denomDF.out,mmmgee=xx)
	out<-list(test=test.out,hypothesis=list(contr.matrix=LL.out,rhs=rr),conf.int=conf.int.out,denomDF=denomDF.out,mmmgee=xx)
	class(out)<-"mmmgeetest"
	out
}



# @title Print Results From mmmgeetest Objects
# @export
print.mmmgeetest<-function(x,...) {
	cat(paste("\n\t","Hypothesis tests for linear contrasts in multiple marginal GEE models\n\n"))
	cat(paste("Statistic:",ifelse(x$test$type=="maximum","Maximum-type","Quadratic form"),x$test$statistic,"statistic\n"))
	cat(paste("Approximation:",c("Multivariate normal","Multivariate t","Chi-squared","F","Scaled F")[x$test$approximation==c("normal","t","Chisq","F","scaled.F")],"\n"))
	cat(paste("Alternative:",c("Undirected","Less","Greater")[x$test$alternative==c("undirected","less","greater")],"\n\nGlobal test:\n"))

	if(x$test$statistic=="Wald" & x$test$type=="quadratic") {
		if(x$test$approximation=="Chisq") cat(paste("Chi-squared = ",format(x$test$global[,1],digits=4),", df = ",x$test$global$df,", p-value = ",format(x$test$global$p,digits=4),sep=""))
		if(x$test$approximation=="F") cat(paste("F = ",format(x$test$global[,1],digits=4),", df1 = ",x$test$global$df1,", df2 = ",round(x$test$global$df2,1),", p-value = ",format(x$test$global$p,digits=4),sep=""))
		if(x$test$approximation=="scaled.F") cat(paste("Scaled F = ",format(x$test$global[,1],digits=4),", df1 = ",x$test$global$df1,", df2-df1+1 = ",round(x$test$global$df2F,1),", p-value = ",format(x$test$global$p,digits=4),sep=""))
	}	
	if(x$test$statistic=="Wald" & x$test$type=="maximum") {
		if(x$test$approximation=="normal") cat(paste(names(x$test$global)[1]," = ",format(x$test$global[,1],digits=4),", p-value = ",format(x$test$global$p,digits=4),sep=""))
		if(x$test$approximation=="t") cat(paste(names(x$test$global)[1]," = ",format(x$test$global[,1],digits=4),", df = ",x$test$global$df2,", p-value = ",format(x$test$global$p,digits=4),sep=""))
	}
	if(x$test$statistic=="Score" & x$test$type=="quadratic") {		
		cat(paste("Chi-squared = ",format(x$test$global[,1],digits=4),", df = ",x$test$global$df,", p-value = ",format(x$test$global$p,digits=4),sep=""))
	}
	if(x$test$statistic=="Score" & x$test$type=="maximum") {		
		cat(paste(names(x$test$global)[1]," = ",format(x$test$global[,1],digits=4),", df = ",x$test$global$df,", p-value = ",format(x$test$global$p,digits=4),sep=""))
	}

	cat("\n\n")
	#if(x$test$statistic=="Wald" & x$test$approximation%in%c("t","F","scaled.F")) {
	#	#cat("\n\n")
	#	if(x$denomDF$denomDFtype=="minnp") dftype<-"calculated as min(n-p)."
	#	if(x$denomDF$denomDFtype=="user") dftype<-"user specified."
	#	cat(paste("Denominator degrees of freedom",dftype))
	#	cat("\n\n")
	#}
	if(!is.null(x$test$closed.test)) {
		cat("Closed testing procedure:\n")
		print(format(x$test$closed.test,digits=4))
		cat("\n")
	}

	if(!is.null(x$conf.int)) {
		if(x$conf.int$alternative=="undirected") sided<-"two-sided" else sided<-"one-sided"
		cat(paste( x$conf.int$conf.level*100,"% simultaneous",sided,"confidence intervals using multivariate",ifelse(x$conf.int$approximation=="normal","normal distribution:","t-distribution:"),"\n"))
		print(format(x$conf.int$conf.int,digits=4))
		cat("\n")
		cat("Note: The confidence intervals correspond to a single-step maximum-type Wald test.\n")
	}
}


#Hilfsfunktionen fuer den closed test
dec2bin<-function(x,digits=8) {
	b<-rep(0,digits)
	for(k in 1:digits) {
		b[k]<-x%%2	
		x<-(x - b[k])/2
	}
	b
}
namen.fn<-function(anz.dim) {
	#create binary numbers from 0 to 2^anz.dim - 1
	x<-0:(2^anz.dim - 1)
	M<-t(sapply(x,dec2bin,digits=anz.dim))
	if(anz.dim>1) out<-M else out<-t(M)  #because in this case the t() in sapply is not required, no sorting
	out 
}

