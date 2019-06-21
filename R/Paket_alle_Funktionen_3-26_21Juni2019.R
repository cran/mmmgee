#The code contains commmented out lines that may in the future be the basis for the implementation of residual weights as alternative to the implemented cluster weights.

#' @title Fit Generalized Estimating Equation Models
#'
#' @description \code{geem2} is a modified version of \code{\link[geeM]{geem}} to fit generalized estimating equation models and to provide objects that can be used
#'	for simultaneous inference across multiple marginal models using \code{\link{mmmgee}} and \code{\link{mmmgee.test}}.
#'	Like \code{geem}, \code{geem2} estimates coefficients and nuisance parameters using generalized estimating equations. The link and variance functions can be
#'	specified by the user and the syntax is similar to \code{\link[stats]{glm}}. 
#'
#' @param formula a formula expression similar to that for \code{\link[stats]{glm}}, of the form \code{response~predictors}. An offset is allowed, as in glm.
#' @param id a vector identifying the clusters. By default, data are assumed to be sorted such that observations in a cluster are in consecutive rows and higher
#'	numbered rows in a cluster are assumed to be later. If \code{NULL}, then each observation is assigned its own cluster.
#' @param waves a non-negative integer vector identifying components of a cluster. For example, this could be a time ordering. If integers are skipped within a cluster,
#'	then dummy rows with weight 0 are added in an attempt to preserve the correlation structure (except if \code{corstr = "exchangeable"} or \code{"independent"}).
#'	This can be skipped by setting \code{nodummy=TRUE}. When assessing missing values, waves are assumed to start at 1, starting at a larger integer is therefore computationally
#'	inefficient.
#' @param data an optional data frame containing the variables in the model.
#' @param family will determine the link and variance functions. The argument can be one of three options: a family object, a character string, or a list of functions. 
#'	For more information on how to use family objects, see \code{\link[stats]{family}}. If the supplied argument is a character string, then the string should
#'	correspond to one of the family objects.
#'	In order to define a link function, a list must be created with the components \code{(LinkFun, VarFun, InvLink, InvLinkDeriv)}, all of which are vectorized functions.
#'	If the components in the list are not named as \code{(LinkFun, VarFun, InvLink, InvLinkDeriv)}, then \code{geem2} assumes that the functions are given in that order.
#'	\code{LinkFun} and \code{VarFun} are the link and variance functions. \code{InvLink} and \code{InvLinkDeriv} are the inverse of the link function and the derivative
#'	of the inverse of the link function and so are decided by the choice of the link function.
#' @param corstr a character string specifying the correlation structure. Allowed structures are: \code{"independence", "exchangeable", "ar1", "m-dependent", "unstructured", "fixed"},
#'	and \code{"userdefined"}. Any unique substring may be supplied. If \code{"fixed"} or \code{"userdefined"}, then \code{corr.mat} must be specified. If \code{"m-dependent"},
#'	then \code{Mv} is relevant.
#' @param Mv for \code{"m-dependent"}, the value for \code{m}.
#' @param weights A vector of weights for the inverse of the scale factor each observation. 
#'	If an observation is assigned weight 0, it is excluded from the calculations of any parameters. 
#'	Observations with a \code{NA} in any variable will be assigned a weight of 0.
#'	Note that weights are defined differently in \code{geem2} and \code{geem}, see details. 
# @param residual.weights A vector of weights for the residual of each observation. This type of weights can be used to account for missing at random data in terms of inverse probability weighting.
#	To do so, the residual weights should be the inverse of the probability that the observation is non-missing, see details. Both, Residual weights and scale weights can be specified in a model.
#' @param corr.mat the correlation matrix for \code{"fixed"}. Matrix should be symmetric with dimensions >= the maximum cluster size.
#'	If the correlation structure is \code{"userdefined"}, then this is a matrix describing which correlations are the same. In that case, all entries have to be integers,
#'	and values less or equal zero indicate a correlation of zero. The information regarding the user-defined structure are extracted from the upper triangle of the provided matrix.
#' @param init.beta an optional vector with the initial values of beta. If not specified, then the intercept will be set to \code{InvLink(mean(response))}.
#'	\code{init.beta} must be specified if not using an intercept.
#' @param init.alpha an optional scalar or vector giving the initial values for the correlation. If provided along with \code{Mv>1} or unstructured correlation,
#'	then the user must ensure that the vector is of the appropriate length.
#' @param init.phi an optional initial scale parameter. If not supplied, initialized to 1.
#' @param scale.fix if set to \code{TRUE}, then the scale parameter is fixed at the value of init.phi. See details.
#' @param nodummy if set to \code{TRUE}, then dummy rows will not be added based on the values in waves.
#' @param sandwich if \code{TRUE}, calculate robust variance.
#' @param useP if set to \code{FALSE}, do not use the n-p correction for dispersion and correlation estimates, as in Liang and Zeger. This can be useful when the number
#'	of observations is small, as subtracting p may yield correlations greater than 1.
#' @param maxit maximum number of iterations.
#' @param tol tolerance in calculation of coefficients.
#' @param restriction either a contrast matrix or a list of a contrast matrix and a right hand side vector, defining a restriction on the regression coefficients. See details.
#' @param conv.criterion convergence criterion, either \code{"ratio"} or \code{"difference"}. The default is \code{"ratio"}, using the relative change in regression coefficient
#'	estimates as convergence criterion, like in \code{geem}.
#'	With \code{"difference"} the maximum absolute difference in regression coefficient estimates is used. The latter is required if some coefficient is 0,
#'	e.g. by estimation under some restriction.
#'
#' @details The function is a modification of \code{\link[geeM]{geem}} from the geeM package, such that
#' 	additional output is returned that is required for the calculation of covariance
#' 	matrix across multiple marginal models. In particular the contributions of each subject to the estimating equation are made available in the output. Internal functions regarding
#'	the calculation of matrix inverses were modified to improve the handling of missing data.
#'
#'	In \code{geem2}, weights are defined as scale weights, similar to most othe software. Note that, in contrast, the current version of \code{geem} (version 0.10.1) uses 
#'	residual weights.
#'
#'	The scale parameter \code{phi} is used in estimating the residual working correlation parameters and in estimating the model based (naiv) covariance matrix of the regression coefficients.
#'	Similar as in most other software, requesting \code{scale.fix=TRUE} only has an impact on the latter, while the working correlation is still estimated using an empirical 
#'	scale factor for the residuals. In contrast, \code{geem} uses the fixed scale factor also when estimating the working correlation.
#'
#'	\code{geem2} allows for estimation of regression coefficients under linear restrictions \eqn{L\beta=r}, where \eqn{L} is a contrast matrix, \eqn{\beta}
#'	the vector of regression coefficients and \eqn{r} a real valued right hand side vector. Using the argument \code{restriction}, \eqn{L} and \eqn{r} can be specified.
#'	If only \eqn{L} is specified, \eqn{r} is assumed as null vector. The functionality is in particular required to calculate the generalized score test for linear hypotheses about
#'	\eqn{\beta}. Use \code{conv.criterion="difference"} if any regression coefficient is restricted to 0.
#'
#' @note The option to fit a model with linear restrictions concerning the coefficients is implemented to enable the calculation of a generalized score test.
#'	It may also be used to obtain estimates of the coefficients under restrictions. The model based and robust variance estimates of the restricted 
#'	coefficient estimates are found in the slots \code{restricted.naiv.var} and \code{restricted.var}, respectively. Note that the variance of estimates restricted to a single value
#'	is supposed to be zero, however the calculated variance estimate may deviate from zero within machine accuracy.
#'
#' @return A list with class \code{geem2}, similar to the output
#' 	of \code{geem} from the geeM package. The additional slot \code{sandwich.args} contains components to calculate the sandwich variance estimator for the
#'	fitted model and across models if
#'	applied in the multiple marginal model framework.
#'
#' @author The \code{geem} function was written by Lee McDaniel and Nick Henderson, modifications for \code{geem2} are by Robin Ristl \email{robin.ristl@@meduniwien.ac.at}
#' @references Lee S. McDaniel, Nicholas C. Henderson, Paul J. Rathouz. Fast pure R implementation of GEE: application of the matrix package. The R journal 5.1 (2013): 181.
#' @seealso \code{\link{mmmgee}}, \code{\link[geeM]{geem}}, \code{\link{mmmgee.test}}
#'
#' @examples
#' data(keratosis)
#' m1<-geem2(clearance~trt,id=id,data=keratosis,family=binomial,corstr="independence")
#' summary(m1)
#' m2<-geem2(pain~trt,id=id,data=keratosis[keratosis$lesion==1,],family=gaussian,corstr="independence")
#' summary(m2)
#' geem2(pain~trt,id=id,data=keratosis[keratosis$lesion==1,],family=gaussian,corstr="exchangeable")
#' #
#' data(datasim)
#' mod1<-geem2(Y.lin~gr.lang+x1,id=id,data=datasim,family="gaussian",corstr="exchangeable")
#' summary(mod1)
#' mod2<-geem2(Y.poi~gr.lang+x2,id=id,data=datasim,family="poisson",corstr="unstructured")
#' summary(mod2)
#' mod3<-geem2(Y.bin~gr.lang+x3,id=id,data=datasim,family="binomial",corstr="user",
#'		corr.mat=matrix(c(1,2,3,0, 2,1,2,3, 3,2,1,2, 0,3,2,1),4,4))
#' summary(mod3)
#'
#' @export
geem2 <- function (formula, id, waves = NULL, data = parent.frame(), family = gaussian, 
    corstr = "independence", Mv = 1, weights = NULL, corr.mat = NULL, 
    init.beta = NULL, init.alpha = NULL, init.phi = 1, scale.fix = FALSE, 
    nodummy = FALSE, sandwich = TRUE, useP = TRUE, maxit = 20, 
    tol = 1e-05, restriction=NULL, conv.criterion=c("ratio","difference")) 
{
    call <- match.call()
    famret <- getfam(family)
    if (inherits(famret, "family")) {
        LinkFun <- famret$linkfun
        InvLink <- famret$linkinv
        VarFun <- famret$variance
        InvLinkDeriv <- famret$mu.eta
    } else {
        LinkFun <- famret$LinkFun
        VarFun <- famret$VarFun
        InvLink <- famret$InvLink
        InvLinkDeriv <- famret$InvLinkDeriv
    }

	#new: allow formula to be passed as character
	#nothing is changed if formula already is a formula object.
	formula<-as.formula(formula)

    #they way scale.fix works has been changed: In geem, the fixed scale was also used as denominator when estimating the working correlation
    #this may, however result in a singular covariance matrix if the fixed scale is smaller than the true scale.
    #Other software packages use the fixed scale only in the calculation of the model based covariance matrix of beta, and 
    #geem2 now does so too.
    #if (scale.fix & is.null(init.phi)) {
    #    stop("If scale.fix=TRUE, then init.phi must be supplied")
    #}
	#NEW:
	if(scale.fix) {
		if(is.null(init.phi)) stop("If scale.fix=TRUE, then init.phi must be supplied")
		phi.fix<-init.phi #this is used in the very end to scale the model based covariance matrix of beta. Otherwise we do not need the fixed scale.
	}
    
    #Modification: Set convergence criterion
    conv.criterion<-match.arg(conv.criterion,choices=c("ratio","difference"))

    #Modification: Check if restrictions are feasible
    #check if restriction is a solvable system, i.e. rank(L)==rank(L,r)
    if(!is.null(restriction)) {
        if(!is.matrix(restriction) & !is.list(restriction)) stop("restriction must be a matrix or a list containing a matrix and a vector.")
        if(is.matrix(restriction)) restriction<-list(L=restriction,r=rep(0,dim(restriction)[1]))
        if(is.list(restriction) & length(restriction)==1) restriction[[2]]<-rep(0,dim(restriction[[1]])[1])
	if(is.list(restriction) & (!is.matrix(restriction[[1]]) | !is.vector(restriction[[2]])) ) stop("restriction must be a matrix or a list containing a matrix and a vector.")
        rankL<-qr(restriction[[1]])$rank
        nrowL<-dim(restriction[[1]])[1]
	if(!is.vector(restriction[[2]]) | length(restriction[[2]])!=nrowL) stop("No vector or incorrect dimension of restriction right hand.")	
	
	#Check if L is of full rank, otherwise the restricted estimation does not work. 
	#Give a warning if the original system was contradictory:
	if(rankL!=qr(cbind(restriction[[1]],restriction[[2]]))$rank) stop("The entered restriction has no solution.")
	#For test like the Tukey test the contrast matrix in the null hypothesis has less than full rank.
	#Still the restriction matrix in the restricted estimation needs to be of full rank.
        if(rankL!=nrowL) {
		L<-restriction[[1]]
	        r<-restriction[[2]]
		keepline<-rep(FALSE,nrowL)
		index<-1:nrowL
		for(i in index) {
			set<-index[keepline]
			keepline[i]<-qr(L[c(set,i),])$rank==(sum(keepline)+1)
		}
		restriction[[1]]<-L[keepline,,drop=FALSE]
		restriction[[2]]<-r[keepline]
		###Warnung nicht ausgeben, stoert nur, wenn Tukey Type tests als Score test berechnet werden.
		###warning("Restriction matrix is not of full rank. A reduced matrix with full rank is used in the restricted estimation.")
	}

        
    }

    useP <- as.numeric(useP)
    dat <- model.frame(formula, data, na.action = na.pass)
    nn <- dim(dat)[1]
    if (typeof(data) == "environment") {
        id <- id
        weights <- weights
        if (is.null(call$weights)) weights <- rep(1, nn)
		#Also here: residualWeightsInUse aktivieren
	  #if (is.null(call$residual.weights)) {
	  #	residual.weights <- rep(1, nn)
	  #	residualWeightsInUse<-FALSE
	  #} else {
	  #	residualWeightsInUse<-TRUE
	  #}
        waves <- waves
    } else {
        if (length(call$id) == 1) {
            subj.col <- which(colnames(data) == call$id)
            if (length(subj.col) > 0) {
                id <- data[, subj.col]
            }
            else {
                id <- eval(call$id, envir = parent.frame())
            }
        } else if (is.null(call$id)) {
            id <- 1:nn
        }
        if (length(call$weights) == 1) {
            weights.col <- which(colnames(data) == call$weights)
            if (length(weights.col) > 0) {
                weights <- data[, weights.col]
            } else {
                weights <- eval(call$weights, envir = parent.frame())
		#this is to handle the case where the argument is weights=x and x=NULL.
		#Because then length(call$weights==1), and weights is evaluated as NULL, but in
		#this case all weights should be 1.
		if (is.null(weights)) weights <- rep(1, nn)
            }
        }
        else if (is.null(call$weights)) {
            weights <- rep.int(1, nn)
        } #else {
		#required if scale.weights are passed as vector in the function call and because the name
		#in the call is scale.weights while in the further function it is weights
		#weights<-scale.weights
        #}
      #New block to allow for residual weights, not in use
	#if (length(call$residual.weights) == 1) {
      #      resweights.col <- which(colnames(data) == call$residual.weights)
      #      if (length(resweights.col) > 0) {
      #          residual.weights <- data[, resweights.col]
      #      } else {
      #          residual.weights <- eval(call$residual.weights, envir = parent.frame())
      #      }
	#	#NEW: residualWeightsInUse logical is used in the generalInverse Function, 
	#	#because with residual weights the inverse is taken from a hypothetical correlation matrix subsuming all observations (also missing ones)
	#	#whereas otherwise the matrix is first reduced to the dimension of the actually observed data and then inverted.
	#	residualWeightsInUse<-TRUE
      #  }
      #  else if (is.null(call$residual.weights)) {
      #      residual.weights <- rep.int(1, nn)
	#	residualWeightsInUse<-FALSE
      #}
        #
        if (length(call$waves) == 1) {
            waves.col <- which(colnames(data) == call$waves)
            if (length(waves.col) > 0) {
                waves <- data[, waves.col]
            } else {
                waves <- eval(call$waves, envir = parent.frame())
            }
        } else if (is.null(call$waves)) {
            waves <- NULL
        }
    }

	#new: Check if there are no NAs in the id vector.
	if(any(is.na(id))) stop("id must not contain NAs.")
	#new: Check length of the id vector. (if id=xx is passed with xx in the parent environment,
	#xx may be of wrong length.
	if(length(id)!=dim(dat)[1]) stop("length of id vector does must match the dimension of the data.")

	#new: id may in principle be a factor, but there is a problem if not all factor levels are used, which may happen
	#if entire clusters are dropped. Therefore we transform it to numeric if it is a factor.
	if(is.factor(id)) id<-as.numeric(id)
    dat$id <- id
    dat$weights <- weights
    #Modification not in use: include residual weights
    #dat$residual.weights <- residual.weights
    #Modification not done: Recode waves such that they always start with 1 and are consecutively numbered
    #wavestemp<-c(2,3,5,3,5,6,2,5,6)
    #as.numeric(as.factor(wavestemp))
    #if(!is.null(waves)) waves<-as.numeric(as.factor(waves))
    dat$waves <- waves
    #Now still required. Waves may not be labelled "A","B","C" for instance, but moved to next if statement
    # if (!is.numeric(dat$waves) & !is.null(dat$waves)) 
     #   stop("waves must be either a non-negative integer vector or NULL") #modified to say non-negative
    #Should waves start with 1? It is inefficient, if they start with a larger number, because dummyrows will be added.
    #But if we shift the numbers there may be a conflict with a user defined correlation structure
	#The na.inds is moved to be after the dummyrows filling, becaus in the dummyrows function
	#nas are also filled
	#by moving it, all na filling can be done in one step
   # na.inds <- NULL
   # if (any(is.na(dat))) {
   #     na.inds <- which(is.na(dat), arr.ind = TRUE)
   # }
    if (!is.null(waves)) {
	    #Check requirements for wave:
		if(any(is.na(waves))) stop("waves must not contain NAs.")
		minwave<-ifelse(is.numeric(dat$waves),min(dat$waves),0)
		if (!is.numeric(dat$waves) | minwave<1 | any(as.integer(dat$waves)!=dat$waves) ) stop("waves must be either a non-negative integer vector or NULL")
		#check if waves start with 1 and issue a warning that it may be more efficient to have them start at 1
		#(because otherwise unneccesary dummy lines are added.
		if(minwave>1) warning("waves should start at 1 for full computational efficiency.")
        dat <- dat[order(id, waves), ]
    } else {
        dat <- dat[order(id), ]
    }
    cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", 
        "unstructured", "fixed", "userdefined")
    cor.match <- charmatch(corstr, cor.vec)
    if (is.na(cor.match)) {
        stop("Unsupported correlation structure")
    }
	#moved here from below:
    if (cor.match == 0) {
        stop("Ambiguous Correlation Structure Specification") #moved to charmatch(corstr, cor.vec) above
    }

    if (!is.null(dat$waves)) {
		#some modifications to prepare dummyrows for all waves.
		#e.g. geem would fill waves 2,4 with a dummyrow resulting in waves 2,3,4
		#we need this to be filled entirely e.g. 1,2,3,4 when four waves are present maximally
		#this is done by the new dummyrows function
		
		maxnwaves<-max(dat$waves)
		nwaves <- table(dat$id)
		
        wavespl <- split(dat$waves, dat$id)
        idspl <- split(dat$id, dat$id)
		incomp <- rep(0, length(wavespl))
		incomp[nwaves<maxnwaves]<-1
      #  maxwave <- rep(0, length(wavespl))
      # 
	#	for (i in 1:length(wavespl)) {
	#		if(length(wavespl[[i]]) < nwaves) 
	#	}
       # for (i in 1:length(wavespl)) {
        #    maxwave[i] <- max(wavespl[[i]]) - min(wavespl[[i]]) + 
        #        1
        #    if (maxwave[i] != length(wavespl[[i]])) {
        #        incomp[i] <- 1
        #    }
        #}
        if (!is.element(cor.match, c(1, 3)) & (sum(incomp) > 0) & !nodummy) {
		#dummyrows was modified to fill all missing rows
            dat <- dummyrows(formula, dat, maxnwaves, wavespl, idspl) #argument incomp removed
            #check:
        }
	
      waves <- dat$waves
     
    }
    #NEW: moved entirely outside the if statement. We always need to reassign id and waves, because the data 
	#have been ordered by id above.
    id <- dat$id
    weights <- dat$weights
    #Modification not in use: include residual weights
    #residual.weights <- dat$residual.weights
    na.inds <- NULL
    if (any(is.na(dat))) {
        na.inds <- which(is.na(dat), arr.ind = TRUE)
		#na.inds[,1] are the row indices of missing values
		#na.inds[,2] are the column indices of missing values
    }
    if (!is.null(na.inds)) {
		#check that with residual weights, there are no missing values in X. Not in use.
		#if(residualWeightsInUse) {
		#	 if(any(is.na(dat[,-1]))) stop("When using residual weights, missing values are only allowed in the outcome variable.")
		#}

        weights[unique(na.inds[, 1])] <- 0 #unique(na.inds[, 1]) are all row numbers of rows that have some missing value
        for (i in unique(na.inds)[, 2]) {  
		#modified to include character variables
            if (is.factor(dat[, i]) | is.character(dat[, i]) ) {
                dat[na.inds[, 1], i] <- levels(as.factor(dat[, i]))[1]
            }
            else {
                dat[na.inds[, 1], i] <- median(dat[, i], na.rm = TRUE)
            }
        }
    }
    includedvec <- weights > 0
    #Modification not in use: Set residual weight to 0, if scale weight is 0.
    #This is important because we later use the sum of residual weights as denominator in the estimation of phi and correlation parameters.
    #residual.weights[!includedvec]<-0

    inclsplit <- split(includedvec, id)
    dropid <- NULL
    allobs <- TRUE
    if (any(!includedvec)) {
        allobs <- FALSE
        for (i in 1:length(unique(id))) {
            if (all(!inclsplit[[i]])) {
                dropid <- c(dropid, i)
            }
        }
    }
    if (length(dropid) > 0) {
        dropind <- which(is.element(id, dropid))
        dat <- dat[-dropind, ]
        includedvec <- includedvec[-dropind]
        weights <- weights[-dropind]
	
     	  #Modification not in use: also perform drop for residual weights
        #residual.weights <- residual.weights[-dropind]
        id <- id[-dropind]
	#Modification: inclsplit also needs to be updated after dropping ids, because it is used in other functions later on.
	inclsplit <- split(includedvec, id)

		#sollte hier nicht auch stehen #waves<-waves[-dropind] #??? #Nein, waves kommt danach nicht mehr vor.
    }
    nn <- dim(dat)[1]
    K <- length(unique(id))
    modterms <- terms(formula)
    X <- model.matrix(formula, dat)
    Y <- model.response(dat)
    offset <- model.offset(dat)
    p <- dim(X)[2]
    if (is.null(offset)) {
        off <- rep(0, nn)
    } else {
        off <- offset
    }
    interceptcol <- apply(X == 1, 2, all)
    linkOfMean <- LinkFun(mean(Y[includedvec])) - mean(off)
    if (any(is.infinite(linkOfMean) | is.nan(linkOfMean))) {
        stop("Infinite or NaN in the link of the mean of responses.  Make sure link function makes sense for these data.")
    }
    if (any(is.infinite(VarFun(mean(Y))) | is.nan(VarFun(mean(Y))))) {
        stop("Infinite or NaN in the variance of the mean of responses.  Make sure variance function makes sense for these data.")
    }
    if (is.null(init.beta)) {
        if (any(interceptcol)) {
            init.beta <- rep(0, dim(X)[2])
            init.beta[which(interceptcol)] <- linkOfMean
        } else {
            stop("Must supply an initial beta if not using an intercept.")
        }
    }
    includedlen <- rep(0, K)
    len <- rep(0, K)
    uniqueid <- unique(id)
    tmpwgt <- as.numeric(includedvec)
    idspl <- ifelse(tmpwgt == 0, NA, id)
    includedlen <- as.numeric(summary(split(Y, idspl, drop = TRUE))[, 1])
    len <- as.numeric(summary(split(Y, id, drop = TRUE))[, 1])
    W <- Diagonal(x = weights) #used anywhere?
	#NEW: in almost all functions we have StdErr %*% sqrtW in sereval computations
	#and they are never used separately. 
	#So better define StdErr initially in the main loop as diag(1/var(mu)) %*% sqrtW
    	#not required any more:
    #sqrtW <- sqrt(W)
    included <- Diagonal(x = (as.numeric(weights > 0)))
    #Modification not in use: add residual weights
    #RW <- Diagonal(x = residual.weights)

	#modification not in use:
	 #If we use residual.weights, the scale.weights for residual.weights==0 are not set to 0 but 1, because the observations are implicitly considered
	 #Note that there must be no missing values in the X matrix! See Preisser, J. S., Lohman, K. K., & Rathouz, P. J. (2002). Performance of weighted estimating equations for longitudinal binary data with drop-outs missing at random. Statistics in medicine, 21(20), 3035-3054.
	 #if(residualWeightsInUse) weights[residual.weights==0]<-1
	 #still, included defined above is the vector which indicates the really included observations. In the other functions this may be problemativ since 
	 #included %*% sqrtW %*% RW is used there.
	
    #modification not in use: used in the estimation of phi:
    #sqrtRW <- sqrt(RW)
    
    if (is.null(init.alpha)) {
        alpha.new <- 0.2
        if (cor.match == 4) {
            alpha.new <- 0.2^(1:Mv)
        }
        else if (cor.match == 5) {
            alpha.new <- rep(0.2, sum(1:(max(len) - 1)))
        }
        else if (cor.match == 7) {
            #alpha.new <- rep(0.2, max(unique(as.vector(corr.mat)))) #unique is not required
		alpha.new <- rep(0.2, max(as.vector(corr.mat)))
        }
    }
    else {
        alpha.new <- init.alpha
    }
    if (is.null(init.phi)) {
        phi <- 1
    }
    else {
        phi <- init.phi
    }
    beta <- init.beta
    StdErr <- Diagonal(nn)
    dInvLinkdEta <- Diagonal(nn)
    #Resid <- Diagonal(nn) #this is never used.

    #Prepare working correlation structure
    #NEW: Now common for all correlation structures:
	tmp <- getBlockDiag(len)
	BlockDiag <- tmp$BDiag
      row.vec <- tmp$row.vec
      col.vec <- tmp$col.vec
	struct.vec<-NULL #argument for getAlphaInvGeneral, that is only none-null for user defined (7)

	#case of independence is not required here, because it will be calculated as diag(1/phi) anew in each iteration.
	#I moved the line alpha.new <- "independent"
	#here, out from the loop, because this needs to be defined once only.
    if (cor.match == 1) { #independence
        #R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
	  #R.alpha.inv <- Diagonal(x = rep.int(1/phi, nn))
            alpha.new <- "independent"
		#new:
		biggest.R.alpha<-Diagonal(x = rep.int(1, max(len)))
        #BlockDiag <- getBlockDiag(len)$BDiag #now common line above
    }
    #else if (cor.match == 2) { #ar1
    #    tmp <- buildAlphaInvAR(len)
    #    a1 <- tmp$a1
    #    a2 <- tmp$a2
    #    a3 <- tmp$a3
    #    a4 <- tmp$a4
    #    row.vec <- tmp$row.vec
    #    col.vec <- tmp$col.vec
    #    BlockDiag <- getBlockDiag(len)$BDiag
    #}
    #else if (cor.match == 3) { 
    #    #tmp <- getBlockDiag(len)
    #    #BlockDiag <- tmp$BDiag #done above
    #    n.vec <- vector("numeric", nn)
    #    index <- c(cumsum(len) - len, nn)
    #    for (i in 1:K) {
    #        n.vec[(index[i] + 1):index[i + 1]] <- rep(includedlen[i], 
    #            len[i])
    #    }
    #}
    else if (cor.match == 4) { #m-dependent
        if (Mv >= max(len)) {
            stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
        }
        #tmp <- getBlockDiag(len)
        #BlockDiag <- tmp$BDiag
        #row.vec <- tmp$row.vec
        #col.vec <- tmp$col.vec
        #R.alpha.inv <- NULL #is this line needed?
    }
    else if (cor.match == 5) { #unstructured
        if (max(len^2 - len)/2 > length(len)) {
            stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
        }
        #tmp <- getBlockDiag(len)
        #BlockDiag <- tmp$BDiag
        #row.vec <- tmp$row.vec
        #col.vec <- tmp$col.vec
    }
    else if (cor.match == 6) { #fixed
        corr.mat <- checkFixedMat(corr.mat, len)
	  alpha.new<-as.vector(corr.mat[upper.tri(corr.mat, diag = FALSE)])
		InverseAndCorrmat<-getAlphaInvGeneral(alpha.new, row.vec, col.vec, len, includedvec, inclsplit, cor.match, struct.vec)
		R.alpha.invFixed<-InverseAndCorrmat[[1]]
		biggest.R.alpha<-InverseAndCorrmat[[2]]
	     	#R.alpha.invFixed <- getAlphaInvGeneral(alpha.new, row.vec, col.vec, len, includedvec, inclsplit, cor.match)
		#later in the code:  R.alpha.inv<-R.alpha.invFixed/phi

        #R.alpha.inv <- as(getAlphaInvFixed(corr.mat, len), "symmetricMatrix")/phi
        #BlockDiag <- getBlockDiag(len)$BDiag
    }
    else if (cor.match == 7) { #user defined structure
        corr.mat <- checkUserMat(corr.mat, len)
        tmp1 <- getUserStructure(corr.mat)
        corr.list <- tmp1$corr.list #needed for updateAlpha
        user.row <- tmp1$row.vec
        user.col <- tmp1$col.vec
        struct.vec <- tmp1$struct.vec
        #tmp2 <- getBlockDiag(len)
        #BlockDiag <- tmp2$BDiag
        #row.vec <- tmp2$row.vec
        #col.vec <- tmp2$col.vec
    }
    #else if (cor.match == 0) {
    #    stop("Ambiguous Correlation Structure Specification") #moved to charmatch(corstr, cor.vec) above
    #}
    #else {
    #    stop("Unsupported Correlation Structure") #this is already checked above
    #}

    stop <- FALSE
    converged <- FALSE
    count <- 0
    beta.old <- beta
    unstable <- FALSE
    phi.old <- phi
    while (!stop) {
        count <- count + 1
        eta <- as.vector(X %*% beta) + off
        mu <- InvLink(eta)
	  #Modification: Include the weights directly here! They are only used as StdErr%*%sqrtW anyways and we do not have to repeat this calculation all over.
        diag(StdErr) <- sqrt(weights/VarFun(mu))

        #old:
	  #diag(StdErr) <- sqrt(1/VarFun(mu))

	  #Modification: Calculate standardized residuals once to be used in all functions (updatePhi, updateAlpha)
	  Resid <- StdErr %*% Diagonal(x = Y - mu)

	  ## NEW: behaviour of scale.fix was modified. We now use a scale parameter to calculate the working correlation, similar to other software.
        ## Only in the final calculation of the model based variance the fixed scale is applied.
        ##if (!scale.fix) {
            phi <- updatePhi(Resid, p, included, useP)
		#before: (Y, mu, VarFun, p, StdErr, included, 
            #    includedlen, sqrtW, RW, sqrtRW, useP) #RW, sqrtRW included
        ##}
        phi.new <- phi

        if (cor.match == 1) {
            R.alpha.inv <- Diagonal(x = rep.int(1/phi, nn))
            #alpha.new <- "independent"
        }
        else if (cor.match == 2) {
            alpha.new <- updateAlphaAR(Resid, phi, p, BlockDiag, included, useP)
			#(Y, mu, VarFun, phi, id, 
      	       #   len, StdErr, p, included, includedlen, includedvec, 
            	  #  allobs, sqrtW, RW, BlockDiag, useP) #RW included
	            #R.alpha.inv <- getAlphaInvAR(alpha.new, a1, a2, a3, 
      	       #   a4, row.vec, col.vec, len, includedvec)/phi #added len and includedvec
        }
        else if (cor.match == 3) {
            alpha.new <- updateAlphaEX(Resid, phi, p, BlockDiag, includedlen, useP)
		#(Y, mu, VarFun, phi, id, 
              #  len, StdErr, Resid, p, BlockDiag, included, includedlen, 
               # sqrtW, RW, useP) #RW included
            #R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
        }
        else if (cor.match == 4) {
            if (Mv == 1) {
                alpha.new <- updateAlphaAR(Resid, phi, p, BlockDiag, included, useP)
			#(Y, mu, VarFun, phi, 
                  #id, len, StdErr, p, included, includedlen, 
                  #includedvec, allobs, sqrtW, RW, BlockDiag, useP) #RW included
            }
            else {
                alpha.new <- updateAlphaMDEP(Resid, phi, p, BlockDiag, Mv, included, includedlen, useP)
				#(Y, mu, VarFun, phi, 
	                  #id, len, StdErr, Resid, p, BlockDiag, Mv, included, 
                  	#includedlen, allobs, sqrtW, RW, useP) #RW included
                if (sum(len > Mv) <= p) {
                  unstable <- TRUE
                }
            }
		#moved to end of fitting process:
            #if (any(alpha.new >= 1)) {
            #    stop <- TRUE
            #    warning("some estimated correlation is greater than 1, stopping.")
            #}

            #R.alpha.inv <- getAlphaInvMDEP(alpha.new, len, row.vec, 
            #    col.vec)/phi
        }
        else if (cor.match == 5) {
            alpha.new <- updateAlphaUnstruc(Resid, phi, len, p, BlockDiag, included, includedlen, useP)
			#(Y, mu, VarFun, phi, 
                #id, len, StdErr, Resid, p, BlockDiag, included, 
                #includedlen, allobs, sqrtW, RW, useP) #RW included
		
		#moved to end of fitting process:
            #if (any(alpha.new >= 1)) {
 	      #        stop <- TRUE
              #  warning("some estimated correlation is greater than 1, stopping.")
            #}
           # R.alpha.inv <- getAlphaInvUnstruc(alpha.new, len, 
            #    row.vec, col.vec, includedvec)/phi #includedvec added for test purposes
        }
        else if (cor.match == 6) {
            R.alpha.inv <- R.alpha.invFixed/phi
        }
        else if (cor.match == 7) {
            alpha.new <- updateAlphaUser(Resid, phi, len, p, BlockDiag, included, useP, user.row, user.col, corr.list)
			#(Y, mu, phi, id, len, 
                	#StdErr, Resid, p, BlockDiag, user.row, user.col, 
                	#corr.list, included, includedlen, allobs, sqrtW, RW, 
                	#useP) #RW included
           # R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, 
            #    user.row, user.col, row.vec, col.vec)/phi
        }
       
	  #calculate R.inverse
        if(cor.match %in% c(2,3,4,5,7)) {
		InverseAndCorrmat<-getAlphaInvGeneral(alpha.new, row.vec, col.vec, len, includedvec, inclsplit, cor.match,  struct.vec)
		R.alpha.inv<-InverseAndCorrmat[[1]]/phi
		biggest.R.alpha<-InverseAndCorrmat[[2]]
	  }

        beta.list <- updateBeta(Y, X, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, weights, restriction) 
				     #YY, XX, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, tol, weights, restriction
                            ###Modification to match function using variance weights #restriction added by Robin
        beta <- beta.list$beta
        phi.old <- phi

        #Modification to ensure convergence is identified under restrictions like beta=0
        #if (max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < 
        #    tol) {
        if ( (conv.criterion=="difference" & max(abs(beta - beta.old)) < tol) | 
		(conv.criterion=="ratio" & max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < tol) ) {
            converged <- TRUE
            stop <- TRUE
        }
        if (count >= maxit) {
            stop <- TRUE
        }
        beta.old <- beta
    }
	#biggest.R.alpha is now in output list of the function that calculates the inverse
   # biggest <- which.max(len)[1]
   # index <- sum(len[1:biggest]) - len[biggest]
   # if (K == 1) {
   #     biggest.R.alpha.inv <- R.alpha.inv
   #     if (cor.match == 6) {
   #         biggest.R.alpha <- corr.mat * phi
    #    }
    #    else {
    #        biggest.R.alpha <- solve(R.alpha.inv)
    #    }
    #}
    #else {
    #    biggest.R.alpha.inv <- R.alpha.inv[(index + 1):(index + 
    #        len[biggest]), (index + 1):(index + len[biggest])]
    #    if (cor.match == 6) {
    #        biggest.R.alpha <- corr.mat[1:len[biggest], 1:len[biggest]] * 
     #           phi
     #   }
     #   else {
     #       biggest.R.alpha <- solve(biggest.R.alpha.inv)
     #   }
    #}
    eta <- as.vector(X %*% beta) + off
    #final update of StdErr, as it is updated in getSandwich and so mmmgee has access to the same StdErr.
    mu <- InvLink(eta)
    diag(StdErr) <- sqrt(weights/VarFun(mu))

    if (sandwich) {
        sandvar.list <- getSandwich(Y, X, eta, mu, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, beta.list$hess, StdErr, dInvLinkdEta, BlockDiag)
		###Modification to match function using variance weights
    }
    else {
        sandvar.list <- list()
        sandvar.list$sandvar <- "no sandwich"
    }
    if (!converged) {
        warning("Did not converge")
    }
    if (unstable) {
        warning("Number of subjects with number of observations >= Mv is very small, some correlations are estimated with very low sample size.")
    }   
    if (is.character(alpha.new)) {
        alpha.new <- 0
    }
	#Moved here to end, so fitting is completed even if alpha>1 occurs.
    if (any(alpha.new >= 1)) {
        warning("some estimated correlation is greater than 1.")
    }
    results <- list()
    results$beta <- as.vector(beta)
	#modified to handle scale.fix as described in the notes
	if(scale.fix) results$phi<-phi.fix else results$phi<-phi
	results$scale.fix<-scale.fix
    results$alpha <- alpha.new
    if (cor.match == 6) {
        results$alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat, 
            1) != 0)])
    }
    results$coefnames <- colnames(X)
    results$niter <- count
    results$converged <- converged
    results$naiv.var <- solve(beta.list$hess) #Note: this needs to be modified in case that residual.weights are used,
    #because then the Hessian is not the naive variance
    #When we perform a weighted estimation, where the weighting matrix W is assumed to be the true covariance matrix of
    #residuals AND we include residual weights
    #in a diagonal matrix V for additional IPTW weighting due to missing data, the estiamte for beta is of the form
    #beta.hat = (X'WVX)^-1 X'WVy
    #Naive varaince estimate means, we believe that var(y)= W^-1
    #then var(beta.hat) = (X'WVX)^-1 X'WV W^-1 VWX((X'WVX)^-1)'
    #Note that X'WVX is NOT symmetric and needs to be transposed when appearing at the end of the variance calculation
    #(W and V are symmetric)
    #Note further, when V is the identity matrix, so no residual weights are applied, 
    #var(beta.hat) = (X'WX)^-1 X'W W^-1 WX((X'WX)^-1)' = (X'WX)^-1 X'WX(X'WX)^-1) = (X'WX)^-1 , the usual Hessian
    #not in use:
    #results$invhessian <- solve(beta.list$hess) #this output may be used by getW, if residual weights are enabled.
    #if(all(residual.weights[weights>0]==1)) {
		#results$naiv.var <- results$invhessian
    #} else {
	#	incl<-diag(included)==1
	#	#Bnaiv = X'WV W^-1 VWX, X' is in fact t(X)%*% dInvLinkdEta in glm
	#	W <- crossprod(sqrtW[incl,incl] %*% StdErr[incl,incl], R.alpha.inv[incl,incl] %*% sqrtW[incl,incl] %*% StdErr[incl,incl]) #W is symmetric
	#	#V=RW, i.e. residual weights diagonal matrix
	#	VWX <- RW[incl,incl] %*% W %*% dInvLinkdEta[incl,incl] %*% X[incl,]
	#	Bnaiv <- crossprod(VWX,solve(W,VWX)) #Need to use [incl,incl], because for zero weight entries the inverse would be Inf.
	#
	#	#results$check<-list(invH1=results$invhessian,invH2=solve(t(VWX)%*%X[incl,])) #identical as should be
	#
	#	results$naiv.var<- results$invhessian%*%Bnaiv%*%t(results$invhessian)
		#compare:
		#hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% 
		#       XX, R.alpha.inv %*% sqrtW %*% StdErr %*% RW %*% dInvLinkdEta %*% 
		#          XX)
		#esteq <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% 
		#        XX, R.alpha.inv %*% sqrtW %*% StdErr %*% RW %*% (YY - mu))
    #}
	#NEW: new handling of scale.fix, estimated phi is used in calculating correlations, but in the end the naive variance is scaled by the assumed scale factor. 
	if(scale.fix) results$naiv.var<-results$naiv.var/phi*phi.fix
    results$var <- sandvar.list$sandvar
	#NEW:
	results$sandwich<-sandwich
    results$call <- call
    results$corr <- cor.vec[cor.match]
    results$clusz <- len
    results$FunList <- famret

    ####Modification:
    ###NOTE: Below "dat <- model.frame ...", X are defined anew as in the beginning. But after the initial definition,
    ###the missing values in dat (and so in X) are set to the first factor class or the median of the respective variable.
    ###The weights for these observations are then set to zero in sqrtW, effectively removing the observations
    ###from all calculations. For MMM GEE we need ths modified X matrix, not the one containing NAs.
    ###Therefore we keep this X and include it to our additional output.
	
    ###FURTHER NOTE: In the original geem code it seems it was overlooked to do the same with Y. So, if a value in Y is missing
    ###the slot y in the geem output contains a vector in which missing values are replaced by the median of the non-missing values.
    ###I correct this below by adding the line "Y <- model.response(dat)"

    X.namod<-X #added by Robin
    Y.namod<-Y #added by Robin

	#output waves, the score test function accesses this slot
	results$waves<-dat$waves

    dat <- model.frame(formula, data, na.action = na.pass)
    X <- model.matrix(formula, dat)

    Y <- model.response(dat) #added by Robin

    results$X <- X
    results$offset <- off
    results$eta <- eta
    results$dropped <- dropid
    results$weights <- weights #summery function accesses this slot
    #results$scale.weights <- weights
    #results$residual.weights <- residual.weights
    #
    results$terms <- modterms
    results$y <- Y
    results$biggest.R.alpha <- biggest.R.alpha  #/phi this is no longer needed, because biggest.R.alpha is now directly provided
    results$formula <- formula

	#Additional output required to update the model with restrictions for the score test
	results$settings<-list(Mv=Mv,corr.mat=corr.mat,init.beta=init.beta,init.alpha=init.alpha,
		init.phi=init.phi,scale.fix=scale.fix,
		nodummy=nodummy,useP=useP)



    ####Additional output required for the multiple marginal models sandwich estimation (B is not required but was included in earlier versions).
	#if(sandwich) {
	#	B<-sandvar.list$numsand
	#} else {
	#	B<-"no sandwich"
	#}
    results$sandw.args<-list(Y=Y.namod, X=X.namod, eta=eta, id=id, R.alpha.inv=R.alpha.inv, 
            phi=phi, InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun, hess=beta.list$hess, 
            StdErr=StdErr, dInvLinkdEta=dInvLinkdEta, BlockDiag=BlockDiag, len=len,
            score=beta.list$esteq)
		#B=B,included=included,includedlen=includedlen,len=len,n.vec=n.vec) #added by Robin, redundancies: results$clusz<-len,  results$eta <- eta

    #Modification: The input restriction should be contained in the output, so we can see from the model object, if the model was fit under some restriction.
    results$restriction<-restriction
    #Modification: In case of restricted estimation, include variance estimate of restricted coefficient estimates in output.
	#NOTE: If residual weights are enabled, check if this does work with residual weights.
    if(!is.null(restriction)) {
	A<-results$naiv.var
	L<-restriction[[1]]
	Id<-diag(1,dim(A)[1])
	MM<-Id-A%*%t(L)%*%solve(L%*%A%*%t(L),L)
	#MMa<-A-A%*%t(L)%*%solve(L%*%A%*%t(L))%*%L%*%A
	Vtildenaiv<-MM%*%A%*%t(MM)
	results$restricted.naiv.var<-forceSymmetric(Vtildenaiv, uplo="U")
	#results$testvar<-MMa%*%sandvar.list$numsand%*%MMa
	if(sandwich) {
		Vtilde<-MM%*%results$var%*%t(MM)
		results$restricted.var<-forceSymmetric(Vtilde, uplo="U")
	} else {
		results$restricted.var<-"no sandwich"
	}
    }

	#just for checking, remove later:
	#results$dat<-dat

    class(results) <- "geem2"
    return(results)
}



####
#Functions from geeM Version 0.10,  modifications added to use scale weights and improve handling of weights and of missing values
buildAlphaInvAR <- function (len) 
{
    nn <- sum(len)
    K <- length(len)
    a1 <- a2 <- a3 <- a4 <- vector("numeric", nn)
    index <- c(cumsum(len) - len, nn)
    for (i in 1:K) {
        if (len[i] > 1) {
            a1[(index[i] + 1):index[i + 1]] <- c(rep(-1, times = len[i] - 
                1), 0)
            a2[(index[i] + 1):index[i + 1]] <- c(0, rep(1, times = len[i] - 
                2), 0)
            a3[(index[i] + 1):index[i + 1]] <- c(1, rep(0, times = len[i] - 
                2), 1)
            a4[(index[i] + 1):index[i + 1]] <- c(rep(0, times = len[i]))
        }
        else if (len[i] == 1) {
            a1[(index[i] + 1):index[i + 1]] <- 0
            a2[(index[i] + 1):index[i + 1]] <- 0
            a3[(index[i] + 1):index[i + 1]] <- 0
            a4[(index[i] + 1):index[i + 1]] <- 1
        }
    }
    a1 <- a1[1:(nn - 1)]
    subdiag.col <- 1:(nn - 1)
    subdiag.row <- 2:nn
    row.vec <- c(subdiag.row, (1:nn), subdiag.col)
    col.vec <- c(subdiag.col, (1:nn), subdiag.row)
    return(list(row.vec = row.vec, col.vec = col.vec, a1 = a1, 
        a2 = a2, a3 = a3, a4 = a4))
}


checkFixedMat <- function (corr.mat, len) 
{
    if (is.null(corr.mat)) {
        stop("corr.mat must be specified if using fixed correlation structure")
    }
    if (dim(corr.mat)[1] < max(len)) {
        stop("Dimensions of corr.mat must be at least as large as largest cluster")
    }
    if (!isSymmetric(corr.mat)) {
        stop("corr.mat must be symmetric")
    }
    if (determinant(corr.mat, logarithm = TRUE)$modulus == -Inf) {
        stop("supplied correlation matrix is not invertible.")
    }
    return(corr.mat[1:max(len), 1:max(len)])
}


checkUserMat <- function (corr.mat, len) 
{
    if (is.null(corr.mat)) {
        stop("corr.mat must be specified if using user defined correlation structure")
    }
    if (dim(corr.mat)[1] < max(len)) {
        stop("corr.mat needs to be at least as long as the maximum cluster size.")
    }
    test.vec <- as.vector(corr.mat)

	#modified: better check for integers, and include checking if there are no NAs
    #if (any(abs(test.vec - round(test.vec)) > .Machine$double.eps)) {
    if (any(is.na(test.vec)) | any(test.vec!=as.integer(test.vec))) {
        stop("entries in corr.mat must be integers.")
    }

	#NEW: At least one upper triangle entry in test.vec must be a positive integer for the further functions to work.
	uppertri<-corr.mat[upper.tri(corr.mat)]
	if(!any(uppertri>0)) stop("at least one entry in the upper triangle of corr.mat must be positive.")

	#I think the check below is actually not required, though it is of course computationally inefficient to have non-consecutive integers.
	#Also it will not work if unique(test.vec) has a length different than min.val:max.val
	#And in fact it is not checked here if the values start at 1.
	#And it is not required, since values <=0 will result in zero correlation, which is a useful feature.
    #max.val <- max(test.vec)
    #min.val <- min(test.vec)
    #if (!all(sort(unique(test.vec)) == min.val:max.val)) {
        #stop("entries in corr.mat must be consecutive integers starting at 1.")
    #}
    return(corr.mat[1:max(len), 1:max(len)])
}


coef.geem2 <- function (object, ...) 
{
    coefs <- object$beta
    names(coefs) <- object$coefnames
    return(coefs)
}


#NEW: This enables the default confint function to be applied
#to geem2 objects and calculate unadjusted Wald confidence intervals.
vcov.geem2 <- function (object, ...) {
	if(object$sandwich) {
		V<-as.matrix(object$var) 
	} else {
		V<-as.matrix(object$naiv.var) 
	}
	colnames(V)<-object$coefnames
	rownames(V)<-object$coefnames
	return(V)
}


#Modified to create ALL dummyrows, This is required with unstructured, user defined and fixed working correlation
#incomp also needs to be modified. 
dummyrows <- function (formula, dat, maxnwaves, wavespl, idspl) #incomp, removed
{
	#drl<-list(formula,dat,incomp,maxwave,wavespl,idspl)
	#formula<-drl[[1]]
	#dat<-drl[[2]]
	#incomp<-drl[[3]]
	#maxwave<-drl[[4]]
	#wavespl<-drl[[5]]
	#idspl<-drl[[6]]

    #missing <- missid <- misswave <- rep(0, sum(maxwave))
	#should we also have minwave?
    #index <- 1
  #  for (i in 1:length(wavespl)) {
  #      wa <- wavespl[[i]]
  #      if (incomp[i] == 1 & length(wa)>1) { #bei length(wa)==1 ist incomp[i]==0
	  #if ( length(wa)>1) {

  #      	index <- index + 1
	#	#it only fills entries after the first present wave!
	#	#was ist, wenn length(wa) == 1 ist???
	#     for (j in 2:length(wa)) {
      #	      wdiff <- wa[j] - wa[j - 1] - 1
      #     	if (wdiff > 0) {
	#               missing[index:(index + wdiff - 1)] <- (wa[j - 1] + 1):(wa[j] - 1)
      #     	    missid[index:(index + wdiff - 1)] <- idspl[[i]][1]
	#           }
      #	      index <- index + wdiff + 1
	#	}
	#  } else {
      #     index <- index + length(wa)
	#  }
	#
    #}
    #dat2 <- as.data.frame(matrix(nrow = sum(maxwave), ncol = dim(dat)[2]))
	#NEW:
	datAdd<-as.data.frame(matrix(nrow = maxnwaves*length(wavespl) - dim(dat)[1], ncol = dim(dat)[2]))
	colnames(datAdd) <- colnames(dat)
	dat2 <- rbind(dat,datAdd)

	#dat2 <- as.data.frame(matrix(nrow = maxnwaves*length(wavespl), ncol = dim(dat)[2])) # (*)
	#colnames(dat2) <- colnames(dat) #(*)
		#(*) this did work with offset, because the colname in the model.frame includes the word "offset"
		#whereas the colnames in the data do not. e.g. ~ offset(x) ==> colname "offset(x)" but this will not work if the formula
		#is applied anwew!
	
	#missid<-rep(NA,dim(dat2)[1])
	missing<-rep(FALSE,dim(dat2)[1])
	indexsetbasis<-indexset<-1:maxnwaves
	for(i in 1:length(wavespl)) {
		fillset<-indexset[!(indexsetbasis%in%wavespl[[i]])]
		missing[fillset]<-TRUE
		#missid[fillset]<-idspl[[i]][1] #this does not work if id is a factor
		dat2$id[fillset]<-idspl[[i]][1] #this works if id is a factor
		indexset<-indexset+maxnwaves
	}
	
    
	#NEW:
	#Now dat2 is filled anew with lines from dat.
    dat2[!missing, ] <- dat
    #dat2$id[missing] <- missid[missing]
    dat2$weights[missing] <- 0
    dat2$waves<-1:maxnwaves #[missing == 1] <- missing[missing > 0]
    #Modification not in use: Include residual.weights
    #dat2$residual.weights[missing] <- 0 #otherwise these entries would be NA which is not allowed for further calculations

	#this is not required any more, because NAs are filled now after the dummyrows are added
   # NAcols <- which(!is.element(names(dat2), c("id", "waves", 
    #    "weights", "residual.weights"))) #include residual.weights in this list
    #for (i in NAcols) {
    #    dat2[missing, i] <- median(dat2[, i], na.rm = TRUE) #what about string variables?
    #}
    #retdat <- model.frame(formula, dat2, na.action = na.pass)
    #retdat$id <- dat2$id
    #retdat$weights <- dat2$weights
    #Modification: Include residual weights
    #retdat$residual.weights <- dat2$residual.weights
    #retdat$waves <- dat2$waves
    #return(retdat)
	return(dat2) #We can now just return dat2, because it is a model.frame, this is necessary to extract the offset
			#by model.offset(dat) later on.
}


family.geem2 <- function (object, ...) 
{
    return(object$FunList)
}


fitted.geem2 <- function (object, ...) 
{
    InvLink <- object$FunList[[if (object$FunList$family == "custom") 
        "InvLink"
    else "linkinv"]]
    return(InvLink(object$eta))
}

#The functions getAlphaInvAR, getAlphaInvEX, getAlphaInvFixed, getAlphaInvMDEP, getAlphaInvUnstruc, getAlphaInvUser
#were removed and replaced by the new function getAlphaInvGeneral


getBlockDiag <- function (len, xvec = NULL) 
{
    K <- length(len)
    if (is.null(xvec)) {
        xvec <- rep.int(1, sum(len^2))
    }
    row.vec <- col.vec <- vector("numeric", sum(len^2))
    add.vec <- cumsum(len) - len
    if (K == 1) {
        index <- c(0, sum(len^2))
    }
    else {
        index <- c(0, (cumsum(len^2) - len^2)[2:K], sum(len^2))
    }
    for (i in 1:K) {
        row.vec[(index[i] + 1):(index[i + 1])] <- rep.int((1:len[i]) + 
            add.vec[i], len[i])
        col.vec[(index[i] + 1):(index[i + 1])] <- rep((1:len[i]) + 
            add.vec[i], each = len[i])
    }
    BlockDiag <- sparseMatrix(i = row.vec, j = col.vec, x = xvec)
    if (!is.null(xvec)) {
        testsymm <- abs(sum(skewpart(BlockDiag)))
        if (testsymm != 0) {
            warning("Correlation matrix is not computed to be exactly symmetric. Taking only the symmetric part.")
        }
    }
    return(list(BDiag = symmpart(BlockDiag), row.vec = row.vec, 
        col.vec = col.vec))
}


getfam <- function (family) 
{
    if (is.character(family)) {
        family <- get(family, mode = "function", envir = parent.frame(2))
    }
    if (is.function(family)) {
        family <- family()
        return(family)
    }
    else if (inherits(family, "family")) {
        return(family)
    }
    else if (is.list(family)) {
        if (length(match(names(family), c("LinkFun", "VarFun", 
            "InvLink", "InvLinkDeriv"))) == 4) {
            famname <- "custom"
            LinkFun <- family$LinkFun
            InvLink <- family$InvLink
            VarFun <- family$VarFun
            InvLinkDeriv <- family$InvLinkDeriv
        }
        else {
            famname <- "custom"
            LinkFun <- family[[1]]
            VarFun <- family[[2]]
            InvLink <- family[[3]]
            InvLinkDeriv <- family[[4]]
        }
        FunList <- list(family = famname, LinkFun = LinkFun, 
            VarFun = VarFun, InvLink = InvLink, InvLinkDeriv = InvLinkDeriv)
        return(FunList)
    }
    else {
        stop("problem with family argument: should be string, family object, or list of functions")
    }
}

###Modification: getSandwich is taken from geeM version 0.08 to use the variance weights
getSandwich <- function (YY, XX, eta, mu, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, hessMat, StdErr, dInvLinkdEta, BlockDiag) {
#id removed from argument list, it is not required because BlockDiag is supplied. mu added.
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    #mu and StdErr are now calculated in the main function and provided to this function.
    #mu <- InvLink(eta)
    #diag(StdErr) <- sqrt(weights/VarFun(mu))
    scoreDiag <- Diagonal(x = YY - mu) #May here add modification to include residual weights 
    BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
    numsand <- crossprod(StdErr %*% dInvLinkdEta %*% XX,
		R.alpha.inv %*% StdErr %*% BlockDiag %*% StdErr %*% R.alpha.inv %*% StdErr %*% dInvLinkdEta %*% XX)
	#crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% 
       # XX, R.alpha.inv %*% sqrtW %*% StdErr %*% BlockDiag %*% 
        #StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% 
        #dInvLinkdEta %*% XX)
    sandvar <- t(solve(hessMat, numsand))
    #now here is a crucial step which is wrong in geeM v0.10.0 and v0.10.1
    #this original line is correct only if we assume t(A)==A. This does not hold with residual weights!
    #with all the transposing done, one is unnecessary anyways, even if t(A)==A. 
    #old line in geeM, ok with scale weights, not with residual weights:
    #sandvar <- t(solve(t(hessMat), sandvar))
    #instead here use the generally true (and simpler) form:
    sandvar <- solve(hessMat, sandvar)
    return(list(sandvar = sandvar, numsand = numsand))
}


getUserStructure <- function (corr.mat) 
{
    ml <- dim(corr.mat)[1]
    row.vec <- NULL
    col.vec <- NULL
    for (i in 2:ml) {
        row.vec <- c(row.vec, 1:(i - 1))
        col.vec <- c(col.vec, rep(i, each = i - 1))
    }
    struct.vec <- corr.mat[cbind(row.vec, col.vec)]
    corr.list <- vector("list", max(struct.vec))
    for (i in 1:max(struct.vec)) {
        corr.list[[i]] <- which(struct.vec == i) 
		#note: the numbers on the diagonal result in empty entries in the list, because, row.vec and col.vec refer to the upper triangle only.
		#Values below 1 do not appear in the list at all, because length=max(struct.vec) and the loop starts at 1.
		#Therefore values below 1 can be used to encode a correlation of 0.
    }
    return(list(corr.list = corr.list, row.vec = row.vec, col.vec = col.vec, 
        struct.vec = struct.vec))
}


model.matrix.geem2 <- function (object, ...) 
{
    return(model.matrix(object$formula, data = model.frame(object)))
}


predict.geem2 <- function (object, newdata = NULL, ...) 
{
    coefs <- object$beta
    if (is.null(newdata)) {
        return(as.vector(object$X %*% object$beta))
    }
    else {
        if (dim(newdata)[2] != length(coefs)) {
            warning("New observations must have the same number of rows as coefficients in the model")
        }
        return(as.vector(newdata %*% object$beta))
    }
}


print.geem2 <- function (x, ...) 
{
    coefdf <- signif(data.frame(x$beta), digits = 4)
    rownames(coefdf) <- x$coefnames
    colnames(coefdf) <- ""
    print(x$call)
    cat("\n", "Coefficients:", "\n")
    print(t(coefdf))
    cat("\n Scale Parameter: ", signif(x$phi, digits = 4), "\n")
    cat("\n Correlation Model: ", x$corr)
    if (!is.element(x$corr, c("independence", "ar1", "exchangeable"))) {
        if (dim(x$biggest.R.alpha)[1] > 4) {
            cat("\n Working Correlation[1:4,1:4]: \n")
            print(as.matrix(round(x$biggest.R.alpha[1:4, 1:4], 
                digits = 4)))
        }
        else {
            cat("\n Working Correlation: \n")
            print(as.matrix(round(x$biggest.R.alpha, digits = 4)))
        }
    }
    else {
        cat("\n Estimated Correlation Parameter: ", signif(x$alpha, 
            digits = 4), "\n")
    }
    cat("\n Number of clusters: ", length(x$clusz), "  Maximum cluster size: ", 
        max(x$clusz), "\n")
    cat(" Number of observations with nonzero weight: ", sum(x$weights != 
        0), "\n")
}


#Modified as summary.geem2 now includes coefficients as matrix, similar to summary.lm
print.summary.geem2 <- function (x, ...) 
{
    #Coefs <- matrix(0, nrow = length(x$coefnames), ncol = 5)
    #rownames(Coefs) <- c(x$coefnames)
    #colnames(Coefs) <- c("Estimates", "Model SE", "Robust SE", 
    #    "wald", "p")
    #Coefs[, 1] <- x$beta
    #Coefs[, 2] <- x$se.model
    #Coefs[, 3] <- x$se.robust
    #Coefs[, 4] <- x$wald.test
    #Coefs[, 5] <- x$p
	Coefs<-x$coefficients
    #print(signif(Coefs, digits = 4))
    #Modification: With signif there are trailing zeros at end of rounded numbers
    #e.g. beta=c(12.1234,0.1234) will be displayed as 12.1200, 0.1234.
    #better to use just print, which will display all numbers up to the maximally shown position
    #e.g. 12.1234, 0.1234
    print(Coefs, digits = 4)

    if(x$restriction) {
        cat("\n Restricted estimates and their standard errors were calculated.\n")
    }
    if (!is.element(x$corr, c("independence", "ar1", "exchangeable"))) {
        if (dim(x$biggest.R.alpha)[1] > 4) {
            cat("\n Working Correlation[1:4,1:4]: \n")
            print(as.matrix(round(x$biggest.R.alpha[1:4, 1:4], 
                digits = 4)))
        }
        else {
            cat("\n Working Correlation: \n")
            print(as.matrix(round(x$biggest.R.alpha, digits = 4)))
        }
    }
    else {
        cat("\n Estimated Correlation Parameter: ", signif(x$alpha, 
            digits = 4), "\n")
    }
    cat(" Correlation Structure: ", x$corr, "\n")
    #modified to include scale.fix case:
    if(x$scale.fix) scaletext<-" Fixed Scale Parameter: " else scaletext<-" Est. Scale Parameter: "
	cat(scaletext, signif(x$phi, digits = 4), 
        "\n")
    cat("\n Number of GEE iterations:", x$niter, "\n")
    cat(" Number of Clusters: ", length(x$clusz), "   Maximum Cluster Size: ", 
        max(x$clusz), "\n")
    cat(" Number of observations with nonzero weight: ", sum(x$weights != 
        0), "\n")
}


summary.geem2 <- function (object, ...) 
{
    Coefs <- matrix(NA, nrow = length(object$beta), ncol = 5)
    Coefs[, 1] <- c(object$beta)
    naive <- !object$sandwich #is.character(object$var)
    #Modification:
    restriction <- !is.null(object$restriction)
    if (!naive && any(diag(object$var) < 0)) {
        naive <- TRUE
        warning("Some elements of robust variance estimate < 0.  Reporting model based SE.")
    }

    if(!restriction) {	
        Coefs[, 2] <- sqrt(diag(object$naiv.var))
        if (naive) {
            Coefs[, 3] <- rep(NA, length(object$beta))
        }
        else {
            Coefs[, 3] <- sqrt(diag(object$var))
        }
    } else {
        Coefs[, 2] <- sqrt(diag(object$restricted.naiv.var))
        if (naive) {
            Coefs[, 3] <- rep(NA, length(object$beta))
        }
        else {
            Coefs[, 3] <- sqrt(diag(object$restricted.var))
        }
    }
    if (naive) {
        Coefs[, 4] <- Coefs[, 1]/Coefs[, 2]
    }
    else {
        Coefs[, 4] <- Coefs[, 1]/Coefs[, 3]
    }
    Coefs[, 5] <- round(2 * pnorm(abs(Coefs[, 4]), lower.tail = FALSE), 
        digits = 8)
    colnames(Coefs) <- c("Estimate", "Model SE", "Robust SE", 
        "Wald", "p")
	#modified, include Coefs directly so users can access the coef table as with lm objects
	#include names
	rownames(Coefs) <- c(object$coefnames)
    summ <- list(
		coefficients=Coefs,
		#beta = Coefs[, 1], se.model = Coefs[, 2], se.robust = Coefs[, 3], wald.test = Coefs[, 4], p = Coefs[, 5],
		 alpha = object$alpha, 
        corr = object$corr, phi = object$phi,
	#include scale.fix
	scale.fix=object$scale.fix,
	niter = object$niter, 
        clusz = object$clusz, 
		#coefnames = object$coefnames, 
		weights = object$weights, 
        biggest.R.alpha = object$biggest.R.alpha,restriction=restriction)
    class(summ) <- "summary.geem2"
    return(summ)
}

#Modification not in use: RW included

updateAlphaAR <- function(Resid, phi, p, BlockDiag, included, useP) {
	#(YY, mu, VarFun, phi, id, len, StdErr, p, included, 
    #includedlen, includedvec, allobs, sqrtW, RW, BlockDiag, useP) 

    #Resid <- StdErr %*% sqrtW %*% Diagonal(x = YY - mu)
    #Modification to match version 0.10.1 of geeM: insert included (it was missing compared to the other update functions
    #Further modification not in use: residual weights RW added
    #Resid <- StdErr %*% included %*% sqrtW %*% RW %*% Diagonal(x = YY - mu)
    #denom <- phi * (sum(band(triu(included %*% BlockDiag %*% 
    #    included, k = 1), k1 = 1, k2 = 1)) - useP * p)
    #Modification not in use: denominator must be the sum of all residual weights, not the number of summands nincl.
    #In the special case of all residual weights equal to 1, it is the same of course.
    #The degrees of freedom adjustement is nincl/(nincl-p) in both cases. 
    #denompart1 <- sum(band(triu(RW%*%BlockDiag%*%RW, k = 1), k1 = 1, k2 = 1)
    #denom <- ifelse(useP==TRUE, phi * denompart1 * (nincl - p) / nincl, phi * denompart1)

    nincl<-sum(band(triu(included %*% BlockDiag %*% included, k = 1), k1 = 1, k2 = 1))
    denom <- ifelse(useP==TRUE, phi * (nincl - p), phi * nincl)
    #Resid <- StdErr %*% Diagonal(x = YY - mu)
    num <- sum(band(triu(Resid %*% BlockDiag %*% Resid, k = 1), k1 = 1, k2 = 1))
    alpha <- num/denom
    return(as.numeric(alpha))
}


###Modification: updateAlphaEX is taken from geeM version 0.08 to use the variance weights
#If in the future enabling residual weights, note:
#Under missing data and with IPTW weighting using residual weights,
#actually we need IPT weights for the joint probability of pair of residuals beeing non-missing.
#Clearly, a standard IPTW approach does not provide this information.
#We use the product of probabilites as approximation, i.e. we multiply the RW weighted residuals.
#Note that the denominator needs to account for the weights, which is not done in geeM v0.10.1.
#It is done here.
#Further note: If the unweighted residuals are used for estimation in an otherwise correctly
#specified IPT weighted model, the estimate of the correlation can be severely biased.

updateAlphaEX <- function(Resid, phi, p, BlockDiag, includedlen, useP) {
    #not in use: Add RW
    #Resid <- StdErr %*% included %*% sqrtW %*% RW %*% Diagonal(x = YY - mu)
		#Resid <- StdErr %*% Diagonal(x = YY - mu)
    #denominator taking into account residual weights
    #denompart1 <- sum(triu(RW %*% BlockDiag %*% RW, k = 1))
    nincl <- sum(includedlen*(includedlen-1))/2 #Total number of summands
    #nincl and denompart1 are identical if all RW are 1
    #denom <- ifelse(useP==TRUE, phi * denompart1 * (nincl - p) / nincl, phi * denompart1)
    denom <- ifelse(useP==TRUE, phi * (nincl - p), phi * nincl)

    #Version 0.08:
    #denom <- phi * (sum(triu(included %*% BlockDiag %*% included, 
    #    k = 1)) - useP * p)
    #This corresponds to geeM 0.10.1:
    #denom <- phi*(sum(includedlen*(includedlen-1))/2 - useP * p)
    alpha <- sum(triu(Resid %*% BlockDiag %*% Resid, k = 1))
    alpha.new <- alpha/denom
    return(alpha.new)
}

#Modification: This function was modified more strongly than the others
#to avoid distinguisihing so many cases.
#not in use: Include RW and use sum of RWs as denominator
updateAlphaMDEP <- function(Resid, phi, p, BlockDiag, m, included, includedlen, useP) {
	#(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, 
	#    m, included, includedlen, allobs, sqrtW, RW, useP) 
    #Resid <- StdErr %*% included %*% sqrtW %*% RW %*% Diagonal(x = YY - mu) #RW included
 	#RW %*% BlockDiag %*% RW

    #Resid <- StdErr %*% Diagonal(x = YY - mu)
    alpha.new <- vector("numeric", m)
    #not in use:
	#RWBlockDiag <- RW %*% BlockDiag %*% RW
	InclBlockDiag <- included %*% BlockDiag %*% included

	#this line was missing in some intermediate version:
    BlockDiag <- Resid %*% BlockDiag %*% Resid
	for (i in 1:m) {
		#for some reason band does not work with the class of symmetric sparse Matrics "dsCMatrix".
		#but RWBlockDiag is automatically defined as dsCMatrix if all RW are 1.
		#BlockDiag is of course also symmetric, but may in general not be recognized as such automatically
		#to be sure the band function works, I add as( ,"dgCMatrix") below:
		
		bandmat <- drop0(band( as(BlockDiag,"dgCMatrix"), i, i)) #drop0 was required for using bandmat@i below to find nincl
		#bandmat <- band( as(BlockDiag,"dgCMatrix"), i, i)
		num <- sum(bandmat)

		Inclbandmat <- drop0(band(  as(InclBlockDiag,"dgCMatrix")  , i, i)) #drop0 is required for using bandmat@i below to find nincl
		#RWbandmat <- drop0(band(  as(RWBlockDiag,"dgCMatrix")  , i, i)) 
		#RWbandmat <- band(  as(RWBlockDiag,"dgCMatrix")  , i, i)

		#denompart1 <- sum(RWbandmat)
		#this results in not counting observations with a residual of 0:
		#nincl <- length(bandmat@i)
		#this is better, include all non-zero weight observations (recall that RW <- 0 if weight == 0):
		nincl <- length(Inclbandmat@i)

		#denom <- ifelse(useP & (sum(includedlen > i) > p), phi * denompart1 * (nincl - p) / nincl,phi * denompart1)
		denom <- ifelse(useP & (sum(includedlen > i) > p), phi * (nincl - p), phi * nincl)

		alpha.new[i] <- num/denom
	}
	#Note: bandmat@i contains the row numbers of non-zero elements (row counting starts with 0,
	#but this is not important here.
    #old (geeM version 0.10.1
    #for (i in 1:m) {
    #    if (sum(includedlen > i) > p) {
    #        bandmat <- drop0(band(BlockDiag, i, i))
    #        if (allobs) {
    #            alpha.new[i] <- sum(bandmat)/(phi * (sum(as.numeric(len > 
    #              i) * (len - i)) - useP * p))
    #        }
    #        else {
    #            alpha.new[i] <- sum(bandmat)/(phi * (length(bandmat@i) - 
    #              useP * p))
    #       }
    #    }
    #   else {
    #        bandmat <- drop0(band(BlockDiag, i, i))
    #        if (allobs) {
    #            alpha.new[i] <- sum(bandmat)/(phi * (sum(as.numeric(len > 
    #              i) * (len - i))))
    #        }
    #       else {
    #            alpha.new[i] <- sum(bandmat)/(phi * length(bandmat@i))
    #        }
    #    }
    #}
    return(alpha.new)
}

#not in use: include RW and use weights in denominator
updateAlphaUnstruc <- function(Resid, phi, len, p, BlockDiag, included, includedlen, useP) {
	#(YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, 
    	#included, includedlen, allobs, sqrtW, RW, useP) 
    #Resid <- StdErr %*% included %*% sqrtW %*% RW %*% Diagonal(x = YY - mu) #RW included

	#Resid <- StdErr %*% Diagonal(x = YY - mu)

    ml <- max(len)
    #New: Product of weights with same structure as product of residuals (needs to be before next line in order to 
    #use original BlockDiag Matrix.
    #RWBlockDiag <- RW %*% BlockDiag %*% RW
    InclBlockDiag <- included %*% BlockDiag %*% included
    BlockDiag <- Resid %*% BlockDiag %*% Resid

    alpha.new <- vector("numeric", sum(1:(ml - 1)))
    lalph <- length(alpha.new)
    row.vec <- NULL
    col.vec <- NULL
    for (i in 2:ml) {
        row.vec <- c(row.vec, 1:(i - 1))
        col.vec <- c(col.vec, rep(i, each = i - 1))
    }
    index <- cumsum(len) - len
    if (useP & sum(includedlen == max(len)) <= p) {
        stop("Number of clusters of largest size is less than p. Consider setting useP=FALSE.")
    }
    for (i in 1:lalph) {
        newrow <- index[which(len >= col.vec[i])] + row.vec[i]
        newcol <- index[which(len >= col.vec[i])] + col.vec[i]
        bdtmp <- BlockDiag[cbind(newrow, newcol)]
        #New denominator, accounting for weights not in use
	#denompart1mat<-RWBlockDiag[cbind(newrow, newcol)]
	denomMat<-InclBlockDiag[cbind(newrow, newcol)]
        #denompart1 <- sum(denompart1mat)
        #nincl <- sum(bdtmp != 0) #Total number of summands
	#better:
	nincl <- sum(denomMat) #Total number of summands, in this way zero residuals do not result in wrong denominator
      denom <- ifelse(useP, phi * (nincl - p), phi * nincl)
	
        #replaced code
        #if (allobs) {
        #    denom <- (phi * (length(newrow) - useP * p))
        #}
        #else {
        #    denom <- (phi * (sum(bdtmp != 0) - useP * p))
        #}
        alpha.new[i] <- sum(bdtmp)/denom
    }
    return(alpha.new)
}

#include RW and use weights in denominator
updateAlphaUser <- function(Resid, phi, len, p, BlockDiag, included, useP, row.vec, col.vec, corr.list) {
    #Resid <- StdErr %*% included %*% sqrtW %*% RW %*% Diagonal(x = YY - 
     #   mu) #RW included
    ml <- max(len)
    #New: Product of weights with same structure as product of residuals (needs to be before next line in order to 
    #use original BlockDiag Matrix.
    #RWBlockDiag <- RW %*% BlockDiag %*% RW
    #use included only:
    InclBlockDiag <- included %*% BlockDiag %*% included
    BlockDiag <- Resid %*% BlockDiag %*% Resid
    alpha.new <- vector("numeric", length(corr.list))
    index <- cumsum(len) - len

    for (i in 1:length(alpha.new)) {
        newrow <- NULL
        newcol <- NULL
		#for corr.list[[i]] that is empty, there is a lot of [0] indexing, but in the end the result is integer(0), which simply
		#behaves like NULL and does not extend the newrow and newcol vectors.
        for (j in 1:length(corr.list[[i]])) {
            newrow <- c(newrow, index[which(len >= col.vec[corr.list[[i]]][j])] + 
                row.vec[corr.list[[i]][j]])
            newcol <- c(newcol, index[which(len >= col.vec[corr.list[[i]]][j])] + 
                col.vec[corr.list[[i]][j]])
        }
        bdtmp <- BlockDiag[cbind(newrow, newcol)]
	  
        #Denominator
	  denomMat<-InclBlockDiag[cbind(newrow, newcol)]
	  nincl <- sum(denomMat) #Total number of summands, in this way zero residuals do not result in wrong denominator
	  #New denominator, accounting for weights
	  #denompart1mat<-RWBlockDiag[cbind(newrow, newcol)]
	  #denompart1 <- sum(denompart1mat)
	  
	  #nincl <- sum(denompart1mat!=0)
	  denom <- ifelse(useP, phi * (nincl - p),phi * nincl)

        #denom <- ifelse(useP, phi * denompart1 * (nincl - p) / nincl,phi * denompart1)
 		#replaced code
        #if (allobs) {
        #    denom <- phi * (length(newrow) - useP * p)
        #}
        #else {
        #    denom <- phi * (sum(bdtmp != 0) - useP * p)
        #}
        alpha.new[i] <- sum(bdtmp)/denom
    }
    return(alpha.new)
}


####Modification: The updateBeta function uses the variance weights like in the original geeM Version 0.08, not in use: also include residual weights
updateBeta <- function(YY, XX, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, 
    StdErr, dInvLinkdEta, tol, weights, restriction)  #restriction=list(L=NULL,r=dim(L)[1])
{
    beta.new <- beta
    conv = FALSE
    for (i in 1:10) {
        eta <- as.vector(XX %*% beta.new) + off
        diag(dInvLinkdEta) <- InvLinkDeriv(eta)
        mu <- InvLink(eta)
        diag(StdErr) <- sqrt(weights/VarFun(mu)) #modified to include weights directly
	  hess <- crossprod(StdErr %*% dInvLinkdEta %*% XX, R.alpha.inv %*% StdErr %*% dInvLinkdEta %*% XX)
        esteq <- crossprod(StdErr %*% dInvLinkdEta %*% XX, R.alpha.inv %*% StdErr %*% (YY - mu))
	  #Modifications not in use: include factor RW in hess and in esteq
        #hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% XX, R.alpha.inv %*% sqrtW %*% StdErr %*% RW %*% dInvLinkdEta %*% XX)
        #esteq <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% XX, R.alpha.inv %*% sqrtW %*% StdErr %*% RW %*% (YY - mu))
        #Modification to enable score test for the null hypotheses L%*%beta==r
        if(is.null(restriction[[1]])) {
            update <- solve(hess, esteq)
        } else {
            L<-restriction[[1]]
            r<-restriction[[2]]
            #InvHess<-solve(hess)
            #lambda<-solve(L%*%InvHess%*%t(L),L%*%beta.new - r + L%*%InvHess%*%esteq)
            #effizienter ab mmmgeee 1.14:
		LInvHess<-t(solve(t(hess),t(L)))
            lambda<-solve(LInvHess%*%t(L),L%*%beta.new - r + LInvHess%*%esteq)
            update <- solve(hess, esteq-t(L)%*%lambda)
        }
        #end of modification
        #Note: the following if statement is not included in the newer version of geeM. I remove it, becaues it causes 
	#problems if beta.new is zero, which may happen under restrictions of the score test
        #if (max(abs(update)/beta.new) < 100 * tol) {
        #    break
        #}
        beta.new <- beta.new + as.vector(update)
    }
    return(list(beta = beta.new, hess = hess, esteq=esteq)) #Modification, the score value esteq is added to the output, needed for score test.
}

#not in use: Also include residual weights RW when estimating phi
updatePhi <- function(Resid, p, included, useP) {
	#before: (YY, mu, VarFun, p, StdErr, included, includedlen, sqrtW, RW, sqrtRW, useP) 
    #nn <- sum(includedlen)
	resid <- diag(Resid)  #Resid is a diagonal matrix, we need the diagonal entries as vector.
    #resid <- diag(StdErr %*% included %*% sqrtW %*% sqrtRW %*% Diagonal(x = YY - mu))
    #Note, not in use: for the residuals we need to calcualte res^2 * res.weight, because res.weight is the inverse probability
    #of seeing this residual.
    #phi <- (1/(sum(included) - useP * p)) * crossprod(resid, 
    #    resid)
    #Modification not in use: The denominator must be based on the sum of the residual weights (RW are set to zero in the beginning of geem2 for excluded observations,
    #so when no residual weights are specified sum(RW) equals nincl)
    nincl<-sum(included)

    #denom <- ifelse(useP==TRUE, sum(RW)*(nincl - p)/nincl, sum(RW))
	denom <- ifelse(useP==TRUE, nincl - p, nincl)

    phi <- crossprod(resid,resid) / denom
    return(as.numeric(phi))
}






####



######## MMM GEE Funktionen

getX<-function(x) x$sandw.args$X

getW<-function(x,biascorr=FALSE) {
	diag(x$sandw.arg$dInvLinkdEta) <- x$sandw.arg$InvLinkDeriv(x$sandw.arg$eta)
	mu <- x$sandw.arg$InvLink(x$sandw.arg$eta)
	diag(x$sandw.arg$StdErr) <- sqrt(1/x$sandw.arg$VarFun(mu))
	#W, alle Teile sollten symmetrisch sein
	#W<-x$sandw.arg$dInvLinkdEta%*%x$sandw.args$StdErr %*% x$sandw.args$sqrtW %*% x$sandw.args$R.alpha.inv %*% x$sandw.args$sqrtW %*% x$sandw.args$StdErr
	#W<-crossprod(x$sandw.args$sqrtW %*% x$sandw.args$StdErr %*% x$sandw.arg$dInvLinkdEta,x$sandw.args$R.alpha.inv %*% x$sandw.args$sqrtW %*% x$sandw.args$StdErr)
	#StdErr now contains sqrtW
	W<-crossprod(x$sandw.args$StdErr %*% x$sandw.arg$dInvLinkdEta,x$sandw.args$R.alpha.inv %*% x$sandw.args$StdErr)

	if(biascorr==TRUE) {
		#Mancl and DeRouen bias corrected covariance matrix estimator
		#H Matrix fuer Mancl und DeRouen Methode (Biaskorrektur) berechnen:
		#Daten muessen dafuer nach id sortiert sein!!!
		K<-length(x$sandw.arg$len)
		Xcluster<-vector(mode="list",length=K)
		startstop<-c(0,cumsum(x$sandw.arg$len))
		for(i in 1:K) Xcluster[[i]]<-x$sandw.arg$X[(startstop[i]+1):startstop[i+1],,drop=FALSE] #wenn man eine Zeile herausholt, wird es als sonst Spaltenvektor interpretiert
		Xdiag<-bdiag(Xcluster)

		#Hier muessen wir aufpassen, weil naiv var mit residual weights nicht direkt die inverse Hesse Matrix ist
		hessInv<-x$naiv.var #solve(x$sandw.args$hess) #with scale weights this is fine
            #hessInv<-solve(x$sandw.args$hess)
		#hessInv<-x$invhessian #output in geem2 not included, but this must be used if residual weights are used.
				
		#Bias Korrektur bei restringierter Schaetzung, der slot restriction enthaelt L und r, wenn es eine Restriktion gibt, sonst ist er NULL
		#x$restriction[[1]] ist L, x$restriction[[2]] ist r
		if(!is.null(x$restriction)) {
			#hessInvAdj<-hessInv-hessInv%*%t(x$restriction[[1]])%*%solve(x$restriction[[1]]%*%hessInv%*%t(x$restriction[[1]]))%*%x$restriction[[1]]%*%hessInv
			hessInvAdj<-hessInv-hessInv%*%t(x$restriction[[1]])%*%solve(x$restriction[[1]]%*%hessInv%*%t(x$restriction[[1]]),x$restriction[[1]]%*%hessInv)
			hessKronecker<-Diagonal(n=K,x=1)%x%hessInvAdj
		} else {
			#hessKronecker<-diag(1,K)%x%hessInv
			hessKronecker<-Diagonal(n=K,x=1)%x%hessInv
		}

		#Vinv<-x$sandw.args$StdErr %*% x$sandw.args$sqrtW %*% x$sandw.args$R.alpha.inv %*% x$sandw.args$sqrtW %*% x$sandw.args$StdErr
		#H<-x$sandw.arg$dInvLinkdEta%*%Xdiag%*%hessKronecker%*%t(Xdiag)%*%x$sandw.arg$dInvLinkdEta%*%Vinv

		#einfacher:
		H<-x$sandw.arg$dInvLinkdEta%*%Xdiag%*%hessKronecker%*%t(Xdiag)%*%W
		I_n<-Diagonal(n=length(x$sandw.arg$id),x=1)
		#W<-W%*%solve(I_n - H)
		#gleichwertig, im Test nicht unbedingt schneller
		W<-t(solve(t(I_n - H),t(W)))
	}
	W
}

getSdiag<-function(x) {
	mu <- x$sandw.arg$InvLink(x$sandw.arg$eta)
	scoreDiag <- Diagonal(x = x$sandw.arg$Y - mu)
	#not used: Residuen mit RW gewichten.
	#scoreDiag <- x$sandw.arg$RW%*%Diagonal(x = x$sandw.arg$Y - mu)
	scoreDiag
}

getID<-function(x) x$sandw.args$id

getHess<-function(x) x$sandw.args$hess

getHessInv<-function(x) x$naiv.var #Achtung hier! Change to x$invhessian if residual weights are in use 
#getHessInv<-function(x) solve(x$sandw.args$hess)
#getHessInv<-function(x) x$invhessian




#' @title Covariance Matrix Estimation for Multiple Marginal GEE Models
#' @description Calculate the covariance matrix for a stacked vector of regression coefficients 
#'	from multiple marginal GEE models fitted with \code{\link{geem2}}.
#' @param x a list of \code{geem} objects fitted with \code{geem2}. The \code{geem} objects must be different models calculated with data from the same subjects.
#'	In particular, the parameter \code{id} in the call to \code{geem2} must refer to the same subjects in each model.
#' @param biascorr logical, if \code{TRUE}, a bias corrected covariance matrix is calculate by extending the
#'	method due to Mancl and DeRouen to multiple models. See references.
#'
#' @return  A list with class \code{mmmgee} containing the following components:
#' \describe{
#'	\item{\code{beta}}{The stacked vector of regression coefficient estimates from the models in \code{x}.}
#'	\item{\code{V}}{The estimated covariance matrix of the regression coefficient estimates.}
#'	\item{\code{A}}{The outer component of \eqn{V=ABA}.}
#'	\item{\code{B}}{The inner component of \eqn{V=ABA}.}
#'	\item{\code{biascorr}}{The value of the input argument \code{biascorr} (logical).}
#'	\item{\code{n}}{A vector with the number of clusters in each model in \code{x}.}
#'	\item{\code{p}}{A vector with number of regression coefficients in each model in \code{x}.}
#' }
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @references Lloyd A. Mancl, Timothy A. DeRouen. A covariance estimator for GEE with improved small sample properties. Biometrics, 2001, 57(1):126-134.
#'	
#' @seealso \code{\link{geem2}}, \code{\link{mmmgee.test}}
#' @examples
#' data(keratosis)
#' m1<-geem2(clearance~trt,id=id,data=keratosis,family=binomial,corstr="independence")
#' m2<-geem2(pain~trt,id=id,data=keratosis[keratosis$lesion==1,],family=gaussian,corstr="independence")
#' mmmgee(x=list(m1,m2),biascorr=TRUE)
#'
#' @export
mmmgee<-function(x,biascorr=FALSE) {
	X.stack<-bdiag(lapply(x,getX))
	W.stack<-bdiag(lapply(x,getW,biascorr=biascorr))
	S.stack<-bdiag(lapply(x,getSdiag))

	#M (Blockdiagonalmatrix mit Einsern, die anzeigen, welche Elemente von S zum selben Cluster gehoeren
	#Die outer Funktion ist vielleicht etwas langsam.
	#Es sollte moeglich sein, es auf die selbe Art wie in geeM:::getBlockDiag zu programmieren,
	#aber so ist es ganz sicher und die Daten muessen auch nicht nach den Clustern geordnet sein.

	ID.stack<-unlist(lapply(x,getID))
	BD.stack<-Matrix(outer(ID.stack,ID.stack,FUN="=="),sparse=TRUE)

	#Sandwich, innerer Teil, funktioniert auch, wenn bei einem oder mehreren Modellen eine oder mehrere ids fehlen. 
	#wegen Struktur 	0		0	1 1  0
	#			a (0, a) =	  a	1 1    a 
	Upart<-t(X.stack)%*%W.stack%*%S.stack
	B<-Upart%*%BD.stack%*%t(Upart)
	#B<-t(X.stack)%*%W.stack%*%S.stack%*%BD.stack%*%S.stack%*%t(W.stack)%*%X.stack

	#Sandwich, aeusserer Teil
	#Hess.stack<-bdiag(lapply(x,getHess))
	#A<-solve(Hess.stack)
	A<-bdiag(lapply(x,getHessInv))

	#ganzes Sandwich, t(A) weil mit residual weights A nicht symmetrisch ist.
	V<-A%*%B%*%t(A)
	#Das Ergebnis kann im Bereich von 10^-15 asymmetrisch sein, das sind numerische Ungenauigkeiten.
	#Daher symmetrisch machen: 
	V<-forceSymmetric(V, uplo="U")

	#Koeffizienten
	beta<-unlist(lapply(x,coef))

	###
	
	###
	#report n and p, so that other functions can calculate degrees of freedom as min(n-p)
	getn<-function(x) length(x$clusz)
	getp<-function(x) length(x$beta)
	n<-sapply(x,getn)
	p<-sapply(x,getp)

	out<-list(beta=beta,V=V,A=A,B=B,biascorr=biascorr,n=n,p=p)
	class(out)<-"mmmgee"
	out
}


#Notwendige Funktionen zur Anwendung von multcomp glht
#NOTE: These need to be added to the namespace file!
vcov.mmmgee<-function(object,...) {
	V<-as.matrix(object$V)
	colnames(V)<-names(object$beta)
	rownames(V)<-names(object$beta)
	return(V)
}
coef.mmmgee<-function(object,...) object$beta


getZeror<-function(L) rep(0,dim(L)[1])
getp<-function(m) length(coef(m))

#does no longer use update(), to avoid problems with different environments that occurred when using update()
mmmscore.test<-function(Modelle,L.list,r.list=NULL,maxit=20,tol=10^(-8),type="quadratic",alternative="undirected",biascorr=FALSE,showwarning=TRUE,...) {
	if(is.null(r.list)) r.list<-lapply(L.list,getZeror)
	#Fit models under restriction
	M<-length(Modelle)
	ModelleRestr<-vector(length=M,mode="list")
	##for(i in 1:M) ModelleRestr[[i]]<-update(Modelle[[i]],restriction=list(L=L.list[[i]],r=r.list[[i]]),conv.criterion="difference",tol=tol)
	##for(i in 1:M) ModelleRestr[[i]]<-eval(substitute(update(Modelle[[i]],restriction=list(L=L.list[[i]],r=r.list[[i]]),conv.criterion="difference",tol=tol)),envir=parent.frame(),enclos=environment(Modelle[i]$formula))

	#NEU:
	for(i in 1:M) {
		ModelleRestr[[i]]<-geem2(
			formula=Modelle[[i]]$sandw.args$Y~Modelle[[i]]$sandw.args$X-1+offset(Modelle[[i]]$offset),
			id=Modelle[[i]]$sandw.args$id,
			waves=Modelle[[i]]$waves,
			family=Modelle[[i]]$FunList,
			corstr=Modelle[[i]]$corr,
			Mv=Modelle[[i]]$settings$Mv,
			weights=Modelle[[i]]$weights,
			corr.mat=Modelle[[i]]$settings$corr.mat,
			init.beta=Modelle[[i]]$settings$init.beta,
			init.alpha=Modelle[[i]]$settings$init.alpha,
			init.phi= Modelle[[i]]$settings$init.phi,
			scale.fix=Modelle[[i]]$settings$scale.fix,
			nodummy=Modelle[[i]]$settings$nodummy, 
			sandwich=Modelle[[i]]$sandwich,
			useP=Modelle[[i]]$settings$useP,
			maxit=maxit,
			tol=tol,
			restriction=list(L=L.list[[i]],r=r.list[[i]]),
			conv.criterion="difference")
	}
	#####
			
	MMM<-mmmgee(ModelleRestr,biascorr=biascorr)
	#Note: MMM$A is usually symmetric, but not in case that residual weights are used!
		
	getScore<-function(x) x$sandw.args$score[,1]
	#getScore(xRestr[[1]])
	score<-unlist(lapply(ModelleRestr,getScore))
	L<-bdiag(L.list)
	df<-dim(L)[1]
	
	#Quadratic form Score Test
	#For the robust Score Test Statistic, see e.g. D. Boos 1992, On Generalized Score Tests, The American Statistician: Vol 46, No 4
	#Similar to the quadratic form Wald test, a generalized inverse is used if the final covariance matrix is below full rank.
	if(type=="quadratic") {
		#T<-as.numeric( t(score) %*% MMM$A%*%t(L)%*%solve( L%*%MMM$V%*%t(L) )%*%L%*%MMM$A %*% score )
		#now allows for generalized inverse:
		ranktol<- .Machine$double.eps
		LVtL<-L%*%MMM$V%*%t(L)
		SVD<-svd(LVtL)
		rank<-sum(abs(SVD$d)>ranktol)
		if(rank<dim(L)[1] & showwarning) {
			warning("Contrast covariance matrix is below full rank. The test is calculated using a generalized inverse.")
		}
		#in diag(x,nrow,ncol) we need to specify nrow and ncol, because in the case that x is a scalar, diag(x) would create
		#a x by x identity matrix.
		GInv<-SVD$u%*%diag(x=ifelse(SVD$d<ranktol,0,1/SVD$d),nrow=length(SVD$d),ncol=length(SVD$d))%*%t(SVD$v) #generalized inverse
		df<-rank
		LUA<-L%*%MMM$A %*% score
		#T<-as.numeric( t(score) %*% MMM$A%*%t(L)%*%solve( L%*%MMM$V%*%t(L),L%*%MMM$A %*% score ) ) #Note: Here A was assumed to be symmetric
		T<-as.numeric( t(LUA) %*% GInv %*% LUA )
		p<-1-pchisq(T,df=df) 
		retval<-list(test=data.frame(Chisq=T,df=df,p=p),restricted.coef=MMM$beta,score=score,type=type)	
	}
	
	#Max Score Test
	if(type=="maximum") {
		Vscore<-L%*%MMM$V%*%t(L)
		Cscore<-cov2cor(as.matrix(Vscore))
		z<-   ( L%*%MMM$A %*% score	) / sqrt(diag(Vscore))
		if(alternative=="undirected") {
			maxz<-max(abs(z))
			lower=-rep(maxz,df)
			upper=rep(maxz,df)
		}
		if(alternative=="greater") {
			maxz<-max(z)
			lower=rep(-Inf,df)
			upper=rep(maxz,df)
		}
		if(alternative=="less") {
			maxz<-min(z)
			lower=rep(maxz,df)
			upper=rep(Inf,df)
		}
		p<-1-pmvnorm(lower=lower,upper=upper,sigma=Cscore,...)
		retval<-list(test=data.frame(maxZ=maxz,df=df,p=p),restricted.coef=MMM$beta,score=score,type=type,alternative=alternative)
		names(retval$test)[1]<-c("Max|Z|","MaxZ","MinZ")[alternative==c("undirected","greater","less")]
	}

	class(retval)<-"mmmscoretest"
	retval
}

#Interne Funktion fuer Verwendung aus mmmgee.test
#x mmmgee object
#L Matrix
#r vector
#df numeric
#The function now allows for L%*%x$V%*%t(L) to be below full rank. This can happen if L is not of full rank,
#e.g. if redundant contrasts are specified,
#or if V is not of full rank which may happen, e.g. if two identical models are in x.
#Note that the man mmmgee.test function checks if Lbeta=r has a solution.
#So only L matrices are possible here that can be collapsed to a smaller full-rank matrix.
#We use a generalized inverse based on singular value decomposition via the function svd
mmmchisq.test.simple<-function(x,L,r=NULL,approximation=c("Chisq","F","scaled.F"),df2=NA,showwarning=TRUE) {
	approximation<-match.arg(approximation,choices=c("Chisq","F","scaled.F"))
	ranktol<- .Machine$double.eps
	LVtL<-as.matrix(L%*%x$V%*%t(L))
	SVD<-svd(LVtL)
	rank<-sum(abs(SVD$d)>ranktol)
	if(rank<dim(L)[1] & showwarning) {
		warning("Contrast covariance matrix is below full rank. The test is calculated using a generalized inverse.")
	}
	GInv<-SVD$u%*%diag(x=ifelse(SVD$d<ranktol,0,1/SVD$d),nrow=length(SVD$d),ncol=length(SVD$d))%*%t(SVD$v) #generalized inverse
	df1<-rank
	Lbeta_r<-L%*%x$beta-r
	T<-as.numeric( t(Lbeta_r)%*%GInv%*%(Lbeta_r) )
	#T<-as.numeric( t(L%*%x$beta-r)%*%solve(L%*%x$V%*%t(L),L%*%x$beta-r) )
	#df1<-dim(L)[1]
	if(approximation=="Chisq") {
		p<-1-pchisq(T,df=df1) 
		retval<-list(test=data.frame(Chisq=T,df=df1,p=p),approximation=approximation)
	} 	
	if(approximation=="F") {
		T<-T/df1
		p<-1-pf(T,df1=df1,df2=df2)
		retval<-list(test=data.frame(F=T,df1=df1,df2=df2,p=p),approximation=approximation)
	}
	if(approximation=="scaled.F") {
		T<-T*(df2-df1+1)/df1/df2
		p<-1-pf(T,df1=df1,df2=df2-df1+1)
		retval<-list(test=data.frame(F=T,df1=df1,df2F=df2-df1+1,p=p),approximation=approximation)
	}
	#class(retval)<-"mmmchisqtest"
	retval
}

mmmmax.test.simple<-function(x,L,r,alternative="undirected",approximation="normal",df2=NULL,...) {
	Vcontr<-L%*%x$V%*%t(L)
	z<-(L%*%x$beta-r)/sqrt(diag(Vcontr))
	df<-dim(L)[1]
	if(alternative=="undirected") {
		maxz<-max(abs(z))
		lower=-rep(maxz,df)
		upper=rep(maxz,df)
	}
	if(alternative=="greater") {
		maxz<-max(z)
		lower=rep(-Inf,df)
		upper=rep(maxz,df)
	}
	if(alternative=="less") {
		maxz<-min(z)
		lower=rep(maxz,df)
		upper=rep(Inf,df)
	}
	if(approximation=="normal" | is.null(df2)) {
		p<-1-pmvnorm(lower=lower,upper=upper,sigma=cov2cor(as.matrix(Vcontr)),...)
		retval<-list(test=data.frame(maxZ=maxz,p=p),type="MaxZ",alternative=alternative)
		names(retval$test)[1]<-c("Max|Z|","MaxZ","MinZ")[alternative==c("undirected","greater","less")]
	} else {
		p<-1-pmvt(lower=lower,upper=upper,df=floor(df2),corr=cov2cor(as.matrix(Vcontr)),...)
		retval<-list(test=data.frame(maxT=maxz,df1=df,df2=df2,p=p),type="MaxT",alternative=alternative)
		names(retval$test)[1]<-c("Max|T|","MaxT","MinT")[alternative==c("undirected","greater","less")]
	}
	retval
}





###Example Data 1
#' Simulated Data Set for a Study of Actinic Keratosis Treatments
#'
#' A data set simulated under the planning assumptions for a study comparing four radiation regimens for a photodynamic 
#' treatment of actinic keratosis. Each patient receives each treatment in a different skin patch and each patch contains four lesions.
#' Variables are patient identifier (id), treatment (trt), lesion identifier within a patient (lesion), the binary outcome clearance success (1=success, 0=no success)
#' reported for each lesion and the metric outcome pain (larger values indicating more pain) reported for each skin patch. The aim of the study is to compare
#' the treatments B, C and D to the reference treatment A in terms of both outcomes.
#'
#' @docType data
#'
#' @usage data(keratosis)
#'
#' @format A data frame.
#'
#' @keywords datasets
#'
#' @examples
#' data(keratosis)
#' head(keratosis)
"keratosis"

###Example Data 2
#' Simulated Data Set With Three Endpoints
#'
#' A data set was simulated with repeated observations of a continous outcome variable Y.lin, a count data outcome Y.poi and a binary outcome Y.bin.
#' In the simulation, the mean of an outcome variable depends on the binary grouping variable gr.lang and one of the continuos predictors x1, x2 and x3.
#' Observations of all outcomes in the same subject, indicated by the variable id, are correlated. Also x1, x2 and x3 are correlated within subjects.
#' Data are independent between subjects. The example code shows how to fit a marginal GEE model for each outcome variable and how to test, with different methods,
#' the null hypothesis that the grouping variable has no effect on any of the three endpoints.
#'
#' @docType data
#'
#' @usage data(datasim)
#'
#' @format A data frame.
#'
#' @keywords datasets
#'
#' @examples
#' data(datasim)
#' head(datasim)
#' mod1<-geem2(Y.lin~gr.lang+x1,id=id,data=datasim,family="gaussian",corstr="exchangeable")
#' mod2<-geem2(Y.poi~gr.lang+x2,id=id,data=datasim,family="poisson",corstr="exchangeable")
#' mod3<-geem2(Y.bin~gr.lang+x3,id=id,data=datasim,family="binomial",corstr="exchangeable")
#' L1<-L2<-L3<-matrix(c(0,1,0),nrow=1)
#' mmmgee.test(list(mod1,mod2,mod3),L=list(L1,L2,L3),statistic="Wald",type="maximum",
#'		biascorr=TRUE,asymptotic=FALSE,closed.test=TRUE)
#' \dontrun{
#' mmmgee.test(list(mod1,mod2,mod3),L=list(L1,L2,L3),statistic="score",closed.test=TRUE)
#' mmmgee.test(list(mod1,mod2,mod3),L=list(L1,L2,L3),statistic="score",type="quadratic",
#'		closed.test=TRUE)
#' mmmgee.test(list(mod1,mod2,mod3),L=list(L1,L2,L3),statistic="Wald",type="quadratic",
#'		biascorr=TRUE,asymptotic=FALSE,closed.test=TRUE)
#' mmmgee.test(list(mod1,mod2,mod3),L=list(L1,L2,L3),statistic="Wald",type="quadratic",
#'		scaled=TRUE,biascorr=TRUE,asymptotic=FALSE,closed.test=TRUE)
#' }
"datasim"









######


