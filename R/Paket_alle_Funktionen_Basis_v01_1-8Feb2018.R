#All functions from geeM

#basiert auf
#source("Z:\\MMMGEEPaket\\alte_geeM_Version\\allfunctions_urspruengliche_Gewichte_score_test_1-5.R")
#source("Z:\\MMMGEEPaket\\geeM_Funktionen\\MMM_GEE_Funktion_2-14experimentell.R")


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
#' @param waves an integer vector identifying components of a cluster. For example, this could be a time ordering. If integers are skipped within a cluster,
#'	then dummy rows with weight 0 are added in an attempt to preserve the correlation structure (except if \code{corstr = "exchangeable"} or \code{"independent"}).
#'	This can be skipped by setting \code{nodummy=TRUE}.
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
#' @param weights A vector of weights for each observation. If an observation has weight 0, it is excluded from the calculations of any parameters.
#'	Observations with a \code{NA} anywhere (even in variables not included in the model) will be assigned a weight of 0. 
#'	Note that weights are defined differently in \code{geem2} and \code{geem}, see details. 
#' @param corr.mat the correlation matrix for \code{"fixed"}. Matrix should be symmetric with dimensions >= the maximum cluster size.
#'	If the correlation structure is \code{"userdefined"}, then this is a matrix describing which correlations are the same.
#' @param init.beta an optional vector with the initial values of beta. If not specified, then the intercept will be set to \code{InvLink(mean(response))}.
#'	\code{init.beta} must be specified if not using an intercept.
#' @param init.alpha an optional scalar or vector giving the initial values for the correlation. If provided along with \code{Mv>1} or unstructured correlation,
#'	then the user must ensure that the vector is of the appropriate length.
#' @param init.phi an optional initial overdispersion parameter. If not supplied, initialized to 1.
#' @param scale.fix if set to \code{TRUE}, then the scale parameter is fixed at the value of init.phi.
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
#' 	matrix across multiple marginal models. In particular the contributions of each subject to the estimating equation are made available in the output.
#'
#'	In \code{geem2}, the square root of the weight of an observation is defined as multiplier of the standard error of that observation in the calculation of the 
#'	estimating equation. Note that in contrast in the current version of \code{geem} (version 0.10.0) the diagonal matrix of weights is
#'	used as multiplier of the working correlation matrix. In case that the weights for all observations in one cluster are identical, both definitions are identical and can be 
#'	understood as weighing the cluster-wise contributions to the estimating equation.
#'
#'	\code{geem2} allows for estimation of regression coefficients under linear restrictions \eqn{C\beta=r}, where \code{C} is a contrast matrix, \eqn{\beta}
#'	the vector of regression coefficients and \code{r} a real values right hand side vector. Using the argument \code{restriction}, \code{C} and \code{r} can be specified.
#'	If only \code{C} is specified, \code{r} is assumed as null vector. The functionality is in particular required to calculate the generalized score test for linear hypotheses about
#'	\eqn{\beta}. Use \code{conv.criterion="difference"} if any regression coefficient is restricted to 0.
#'
#' @note The option to fit a model with restrictions concerning the coefficients is implemented to enable the calculation of a generalized score test.
#'	It may also be used to obtain estimates of the coefficients under restrictions, but note that currently the calculation of the covariance matrix is not adapted and the 
#'	covariance matrix and standard errors in the output are not valid estimates of the actual quantities.
#'	(The underlying series expansion is incorrect in the restricted model as the score vector is in general not zero at the location of
#'	the restricted parameter estimates.)
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
#' data(lesions)
#' m1<-geem2(clearance~trt,id=id,data=lesions,family=binomial,corstr="independence")
#' summary(m1)
#' m2<-geem2(pain~trt,id=id,data=lesions[lesions$lesion==1,],family=gaussian,corstr="independence")
#' summary(m2)
#' geem2(pain~trt,id=id,data=lesions[lesions$lesion==1,],family=gaussian,corstr="exchangeable")
#'
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
    }
    else {
        LinkFun <- famret$LinkFun
        VarFun <- famret$VarFun
        InvLink <- famret$InvLink
        InvLinkDeriv <- famret$InvLinkDeriv
    }
    if (scale.fix & is.null(init.phi)) {
        stop("If scale.fix=TRUE, then init.phi must be supplied")
    }

    #Modification: Set convergence criterion
    conv.criterion<-match.arg(conv.criterion,choices=c("ratio","difference"))

    #Modification: Check if restrictions are feasible
    #still to add: check if restriction is a solvable system, i.e. rank(L)==rank(L,r)
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
        if (is.null(call$weights)) 
            weights <- rep(1, nn)
        waves <- waves
    }
    else {
        if (length(call$id) == 1) {
            subj.col <- which(colnames(data) == call$id)
            if (length(subj.col) > 0) {
                id <- data[, subj.col]
            }
            else {
                id <- eval(call$id, envir = parent.frame())
            }
        }
        else if (is.null(call$id)) {
            id <- 1:nn
        }
        if (length(call$weights) == 1) {
            weights.col <- which(colnames(data) == call$weights)
            if (length(weights.col) > 0) {
                weights <- data[, weights.col]
            }
            else {
                weights <- eval(call$weights, envir = parent.frame())
            }
        }
        else if (is.null(call$weights)) {
            weights <- rep.int(1, nn)
        }
        if (length(call$waves) == 1) {
            waves.col <- which(colnames(data) == call$waves)
            if (length(waves.col) > 0) {
                waves <- data[, waves.col]
            }
            else {
                waves <- eval(call$waves, envir = parent.frame())
            }
        }
        else if (is.null(call$waves)) {
            waves <- NULL
        }
    }
    dat$id <- id
    dat$weights <- weights
    dat$waves <- waves
    if (!is.numeric(dat$waves) & !is.null(dat$waves)) 
        stop("waves must be either an integer vector or NULL")
    na.inds <- NULL
    if (any(is.na(dat))) {
        na.inds <- which(is.na(dat), arr.ind = T)
    }
    if (!is.null(waves)) {
        dat <- dat[order(id, waves), ]
    }
    else {
        dat <- dat[order(id), ]
    }
    cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", 
        "unstructured", "fixed", "userdefined")
    cor.match <- charmatch(corstr, cor.vec)
    if (is.na(cor.match)) {
        stop("Unsupported correlation structure")
    }
    if (!is.null(dat$waves)) {
        wavespl <- split(dat$waves, dat$id)
        idspl <- split(dat$id, dat$id)
        maxwave <- rep(0, length(wavespl))
        incomp <- rep(0, length(wavespl))
        for (i in 1:length(wavespl)) {
            maxwave[i] <- max(wavespl[[i]]) - min(wavespl[[i]]) + 
                1
            if (maxwave[i] != length(wavespl[[i]])) {
                incomp[i] <- 1
            }
        }
        if (!is.element(cor.match, c(1, 3)) & (sum(incomp) > 
            0) & !nodummy) {
            dat <- dummyrows(formula, dat, incomp, maxwave, wavespl, 
                idspl)
            id <- dat$id
            waves <- dat$waves
            weights <- dat$weights
        }
    }
    if (!is.null(na.inds)) {
        weights[unique(na.inds[, 1])] <- 0
        for (i in unique(na.inds)[, 2]) {
            if (is.factor(dat[, i])) {
                dat[na.inds[, 1], i] <- levels(dat[, i])[1]
            }
            else {
                dat[na.inds[, 1], i] <- median(dat[, i], na.rm = T)
            }
        }
    }
    includedvec <- weights > 0
    inclsplit <- split(includedvec, id)
    dropid <- NULL
    allobs <- T
    if (any(!includedvec)) {
        allobs <- F
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
        id <- id[-dropind]
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
    }
    else {
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
        }
        else {
            stop("Must supply an initial beta if not using an intercept.")
        }
    }
    includedlen <- rep(0, K)
    len <- rep(0, K)
    uniqueid <- unique(id)
    tmpwgt <- as.numeric(includedvec)
    idspl <- ifelse(tmpwgt == 0, NA, id)
    includedlen <- as.numeric(summary(split(Y, idspl, drop = T))[, 
        1])
    len <- as.numeric(summary(split(Y, id, drop = T))[, 1])
    W <- Diagonal(x = weights)
    sqrtW <- sqrt(W)
    included <- Diagonal(x = (as.numeric(weights > 0)))
    if (is.null(init.alpha)) {
        alpha.new <- 0.2
        if (cor.match == 4) {
            alpha.new <- 0.2^(1:Mv)
        }
        else if (cor.match == 5) {
            alpha.new <- rep(0.2, sum(1:(max(len) - 1)))
        }
        else if (cor.match == 7) {
            alpha.new <- rep(0.2, max(unique(as.vector(corr.mat))))
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
    Resid <- Diagonal(nn)
    if (cor.match == 1) {
        R.alpha.inv <- Diagonal(x = rep.int(1, nn))/phi
        BlockDiag <- getBlockDiag(len)$BDiag
    }
    else if (cor.match == 2) {
        tmp <- buildAlphaInvAR(len)
        a1 <- tmp$a1
        a2 <- tmp$a2
        a3 <- tmp$a3
        a4 <- tmp$a4
        row.vec <- tmp$row.vec
        col.vec <- tmp$col.vec
        BlockDiag <- getBlockDiag(len)$BDiag
    }
    else if (cor.match == 3) {
        tmp <- getBlockDiag(len)
        BlockDiag <- tmp$BDiag
        n.vec <- vector("numeric", nn)
        index <- c(cumsum(len) - len, nn)
        for (i in 1:K) {
            n.vec[(index[i] + 1):index[i + 1]] <- rep(includedlen[i], 
                len[i])
        }
    }
    else if (cor.match == 4) {
        if (Mv >= max(len)) {
            stop("Cannot estimate that many parameters: Mv >=  max(clustersize)")
        }
        tmp <- getBlockDiag(len)
        BlockDiag <- tmp$BDiag
        row.vec <- tmp$row.vec
        col.vec <- tmp$col.vec
        R.alpha.inv <- NULL
    }
    else if (cor.match == 5) {
        if (max(len^2 - len)/2 > length(len)) {
            stop("Cannot estimate that many parameters: not enough subjects for unstructured correlation")
        }
        tmp <- getBlockDiag(len)
        BlockDiag <- tmp$BDiag
        row.vec <- tmp$row.vec
        col.vec <- tmp$col.vec
    }
    else if (cor.match == 6) {
        corr.mat <- checkFixedMat(corr.mat, len)
        R.alpha.inv <- as(getAlphaInvFixed(corr.mat, len), "symmetricMatrix")/phi
        BlockDiag <- getBlockDiag(len)$BDiag
    }
    else if (cor.match == 7) {
        corr.mat <- checkUserMat(corr.mat, len)
        tmp1 <- getUserStructure(corr.mat)
        corr.list <- tmp1$corr.list
        user.row <- tmp1$row.vec
        user.col <- tmp1$col.vec
        struct.vec <- tmp1$struct.vec
        tmp2 <- getBlockDiag(len)
        BlockDiag <- tmp2$BDiag
        row.vec <- tmp2$row.vec
        col.vec <- tmp2$col.vec
    }
    else if (cor.match == 0) {
        stop("Ambiguous Correlation Structure Specification")
    }
    else {
        stop("Unsupported Correlation Structure")
    }
    stop <- F
    converged <- F
    count <- 0
    beta.old <- beta
    unstable <- F
    phi.old <- phi
    while (!stop) {
        count <- count + 1
        eta <- as.vector(X %*% beta) + off
        mu <- InvLink(eta)
        diag(StdErr) <- sqrt(1/VarFun(mu))
        if (!scale.fix) {
            phi <- updatePhi(Y, mu, VarFun, p, StdErr, included, 
                includedlen, sqrtW, useP)
        }
        phi.new <- phi
        if (cor.match == 2) {
            alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, id, 
                len, StdErr, p, included, includedlen, includedvec, 
                allobs, sqrtW, BlockDiag, useP)
            R.alpha.inv <- getAlphaInvAR(alpha.new, a1, a2, a3, 
                a4, row.vec, col.vec)/phi
        }
        else if (cor.match == 3) {
            alpha.new <- updateAlphaEX(Y, mu, VarFun, phi, id, 
                len, StdErr, Resid, p, BlockDiag, included, includedlen, 
                sqrtW, useP)
            R.alpha.inv <- getAlphaInvEX(alpha.new, n.vec, BlockDiag)/phi
        }
        else if (cor.match == 4) {
            if (Mv == 1) {
                alpha.new <- updateAlphaAR(Y, mu, VarFun, phi, 
                  id, len, StdErr, p, included, includedlen, 
                  includedvec, allobs, sqrtW, BlockDiag, useP)
            }
            else {
                alpha.new <- updateAlphaMDEP(Y, mu, VarFun, phi, 
                  id, len, StdErr, Resid, p, BlockDiag, Mv, included, 
                  includedlen, allobs, sqrtW, useP)
                if (sum(len > Mv) <= p) {
                  unstable <- T
                }
            }
            if (any(alpha.new >= 1)) {
                stop <- T
                warning("some estimated correlation is greater than 1, stopping.")
            }
            R.alpha.inv <- getAlphaInvMDEP(alpha.new, len, row.vec, 
                col.vec)/phi
        }
        else if (cor.match == 5) {
            alpha.new <- updateAlphaUnstruc(Y, mu, VarFun, phi, 
                id, len, StdErr, Resid, p, BlockDiag, included, 
                includedlen, allobs, sqrtW, useP)
            if (any(alpha.new >= 1)) {
                stop <- T
                warning("some estimated correlation is greater than 1, stopping.")
            }
            R.alpha.inv <- getAlphaInvUnstruc(alpha.new, len, 
                row.vec, col.vec)/phi
        }
        else if (cor.match == 6) {
            R.alpha.inv <- R.alpha.inv * phi.old/phi
        }
        else if (cor.match == 7) {
            alpha.new <- updateAlphaUser(Y, mu, phi, id, len, 
                StdErr, Resid, p, BlockDiag, user.row, user.col, 
                corr.list, included, includedlen, allobs, sqrtW, 
                useP)
            R.alpha.inv <- getAlphaInvUser(alpha.new, len, struct.vec, 
                user.row, user.col, row.vec, col.vec)/phi
        }
        else if (cor.match == 1) {
            R.alpha.inv <- Diagonal(x = rep.int(1/phi, nn))
            alpha.new <- "independent"
        }
        beta.list <- updateBeta(Y, X, beta, off, InvLinkDeriv, 
            InvLink, VarFun, R.alpha.inv, StdErr, dInvLinkdEta, 
            tol, sqrtW, restriction) #W, included) 
                            ###Modification to match function using variance weights #restriction added by Robin
        beta <- beta.list$beta
        phi.old <- phi

        #Modification to ensure convergence is identified under restrictions like beta=0
        #if (max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < 
        #    tol) {
        if (max(abs((beta - beta.old)/(beta.old + .Machine$double.eps))) < 
            tol | (conv.criterion=="difference" & max(abs(beta - beta.old)) < tol)) {
            converged <- T
            stop <- T
        }
        if (count >= maxit) {
            stop <- T
        }
        beta.old <- beta
    }
    biggest <- which.max(len)[1]
    index <- sum(len[1:biggest]) - len[biggest]
    if (K == 1) {
        biggest.R.alpha.inv <- R.alpha.inv
        if (cor.match == 6) {
            biggest.R.alpha <- corr.mat * phi
        }
        else {
            biggest.R.alpha <- solve(R.alpha.inv)
        }
    }
    else {
        biggest.R.alpha.inv <- R.alpha.inv[(index + 1):(index + 
            len[biggest]), (index + 1):(index + len[biggest])]
        if (cor.match == 6) {
            biggest.R.alpha <- corr.mat[1:len[biggest], 1:len[biggest]] * 
                phi
        }
        else {
            biggest.R.alpha <- solve(biggest.R.alpha.inv)
        }
    }
    eta <- as.vector(X %*% beta) + off
    if (sandwich) {
        sandvar.list <- getSandwich(Y, X, eta, id, R.alpha.inv, 
            phi, InvLinkDeriv, InvLink, VarFun, beta.list$hess, 
            StdErr, dInvLinkdEta, BlockDiag, sqrtW) #W, included) ###Modification to match function using variance weights
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

    ####Modification:
    ###NOTE: Below "dat <- model.frame ...", X are defined anew as in the beginning. But after the initial definition,
    ###the missing values in dat (and so in X) are set to the first factor class or the median of the respective variable.
    ###The weights for these observations are then set to zero in sqrtW, effectively removing the observations
    ###from all calculations. For MMM GEE we need ths modified X matrix, not the one containing NAs.
    ###Therefor we keep this X and include it to our additional output.
	
    ###FURTHER NOTE: In the original geem code it seems it was overlooked to do the same with Y. So, if a value in Y is missing
    ###the slot y in the geem output contains a vector in which missing values are replaced by the median of the non-missing values.
    ###I correct this below by adding the line "Y <- model.response(dat)"

    X.namod<-X #added by Robin
    Y.namod<-Y #added by Robin

    dat <- model.frame(formula, data, na.action = na.pass)
    X <- model.matrix(formula, dat)

    Y <- model.response(dat) #added by Robin

    if (is.character(alpha.new)) {
        alpha.new <- 0
    }
    results <- list()
    results$beta <- as.vector(beta)
    results$phi <- phi
    results$alpha <- alpha.new
    if (cor.match == 6) {
        results$alpha <- as.vector(triu(corr.mat, 1)[which(triu(corr.mat, 
            1) != 0)])
    }
    results$coefnames <- colnames(X)
    results$niter <- count
    results$converged <- converged
    results$naiv.var <- solve(beta.list$hess)
    results$var <- sandvar.list$sandvar
    results$call <- call
    results$corr <- cor.vec[cor.match]
    results$clusz <- len
    results$FunList <- famret
    results$X <- X
    results$offset <- off
    results$eta <- eta
    results$dropped <- dropid
    results$weights <- weights
    results$terms <- modterms
    results$y <- Y
    results$biggest.R.alpha <- biggest.R.alpha/phi
    results$formula <- formula

    ####Additional output required for the multiple marginal models sandwich estimation
    results$sandw.args<-list(Y=Y.namod, X=X.namod, eta=eta, id=id, R.alpha.inv=R.alpha.inv, 
            phi=phi, InvLinkDeriv=InvLinkDeriv, InvLink=InvLink, VarFun=VarFun, hess=beta.list$hess, 
            StdErr=StdErr, dInvLinkdEta=dInvLinkdEta, BlockDiag=BlockDiag, sqrtW=sqrtW, len=len,
            B=sandvar.list$numsand, score=beta.list$esteq)
		#included=included,includedlen=includedlen,len=len,n.vec=n.vec) #added by Robin, redundancies: results$clusz<-len,  results$eta <- eta

    #Modification: The input restriction should be contained in the output, so we can see from the model object, if the model was fit under some restriction.
    results$restriction<-restriction

    class(results) <- "geem2"
    return(results)
}



####
#All functions from geeM Version 0.10
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
    if (determinant(corr.mat, logarithm = T)$modulus == -Inf) {
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
    if (any(abs(test.vec - round(test.vec)) > .Machine$double.eps)) {
        stop("entries in corr.mat must be integers.")
    }
    max.val <- max(test.vec)
    min.val <- min(test.vec)
    if (!all(sort(unique(test.vec)) == min.val:max.val)) {
        stop("entries in corr.mat must be consecutive integers starting at 1.")
    }
    return(corr.mat[1:max(len), 1:max(len)])
}


coef.geem2 <- function (object, ...) 
{
    coefs <- object$beta
    names(coefs) <- object$coefnames
    return(coefs)
}


dummyrows <- function (formula, dat, incomp, maxwave, wavespl, idspl) 
{
    missing <- missid <- misswave <- rep(0, sum(maxwave))
    index <- 1
    for (i in 1:length(wavespl)) {
        wa <- wavespl[[i]]
        index <- index + 1
        for (j in 2:length(wa)) {
            wdiff <- wa[j] - wa[j - 1] - 1
            if (wdiff > 0) {
                missing[index:(index + wdiff - 1)] <- (wa[j - 
                  1] + 1):(wa[j] - 1)
                missid[index:(index + wdiff - 1)] <- idspl[[i]][1]
            }
            index <- index + wdiff + 1
        }
    }
    dat2 <- as.data.frame(matrix(nrow = sum(maxwave), ncol = dim(dat)[2]))
    colnames(dat2) <- colnames(dat)
    dat2[missing == 0, ] <- dat
    dat2$id[missing > 0] <- missid[missing > 0]
    dat2$weights[missing > 0] <- 0
    dat2$waves[missing > 0] <- missing[missing > 0]
    NAcols <- which(!is.element(names(dat2), c("id", "waves", 
        "weights")))
    for (i in NAcols) {
        dat2[missing > 0, i] <- median(dat2[, i], na.rm = TRUE)
    }
    retdat <- model.frame(formula, dat2, na.action = na.pass)
    retdat$id <- dat2$id
    retdat$weights <- dat2$weights
    retdat$waves <- dat2$waves
    return(retdat)
}


family.geem2 <- function (object, ...) 
{
    return(object$FunList)
}

#Diese Funktion ist im geem Paket enthalten, wird aber nie verwendet.
#fillMatList <- function (real.sizes) 
#{
#}


fitted.geem2 <- function (object, ...) 
{
    InvLink <- object$FunList[[if (object$FunList$family == "custom") 
        "InvLink"
    else "linkinv"]]
    return(InvLink(object$eta))
}


getAlphaInvAR <- function (alpha.new, a1, a2, a3, a4, row.vec, col.vec) 
{
    corr.vec <- c(alpha.new * a1/(1 - alpha.new^2), ((1 + alpha.new^2) * 
        a2 + a3)/(1 - alpha.new^2) + a4, alpha.new * a1/(1 - 
        alpha.new^2))
    return(as(sparseMatrix(i = row.vec, j = col.vec, x = corr.vec), 
        "symmetricMatrix"))
}


getAlphaInvEX <- function (alpha.new, diag.vec, BlockDiag) 
{
    return(as(BlockDiag %*% Diagonal(x = (-alpha.new/((1 - alpha.new) * 
        (1 + (diag.vec - 1) * alpha.new)))) + Diagonal(x = ((1 + 
        (diag.vec - 2) * alpha.new)/((1 - alpha.new) * (1 + (diag.vec - 
        1) * alpha.new)) + alpha.new/((1 - alpha.new) * (1 + 
        (diag.vec - 1) * alpha.new)))), "symmetricMatrix"))
}


getAlphaInvFixed <- function (mat, len) 
{
    K <- length(len)
    mat.sizes <- sort(unique(len))
    mat.inverses <- list()
    sl2 <- sum(len^2)
    corr.vec <- vector("numeric", sl2)
    index <- vector("numeric", K + 1)
    index[1] <- 0
    index[2:K] <- (cumsum(len^2) - len^2)[2:K]
    index[K + 1] <- sl2
    for (i in 1:length(mat.sizes)) {
        tmp <- mat[1:mat.sizes[i], 1:mat.sizes[i]]
        mat.inverses[[i]] <- as.vector(solve(tmp))
    }
    mat.finder <- match(len, mat.sizes)
    corr.vec <- unlist(mat.inverses[mat.finder])
    return(as(getBlockDiag(len, corr.vec)$BDiag, "symmetricMatrix"))
}


getAlphaInvMDEP <- function (alpha.new, len, row.vec, col.vec) 
{
    K <- length(len)
    N <- sum(len)
    m <- length(alpha.new)
    mat.sizes <- sort(unique(len))
    corr.vec <- vector("numeric", sum(len^2))
    mat.inverses <- list()
    index <- c(0, (cumsum(len^2) - len^2)[2:K], sum(len^2))
    for (i in 1:length(mat.sizes)) {
        if (mat.sizes[i] == 1) {
            mat.inverses[[i]] <- 1
        }
        else {
            mtmp <- min(m, mat.sizes[i] - 1)
            a1 = list()
            a1[[1]] <- rep(1, mat.sizes[i])
            for (j in 1:mtmp) {
                a1[[j + 1]] <- rep(alpha.new[j], mat.sizes[i] - 
                  j)
            }
            tmp <- bandSparse(mat.sizes[i], k = c(0:mtmp), diagonals = a1, 
                symmetric = T)
            mat.inverses[[i]] <- as.vector(solve(tmp))
        }
    }
    mat.finder <- match(len, mat.sizes)
    corr.vec <- unlist(mat.inverses[mat.finder])
    return(as(sparseMatrix(i = row.vec, j = col.vec, x = corr.vec), 
        "symmetricMatrix"))
}


getAlphaInvUnstruc <- function (alpha.new, len, row.vec, col.vec) 
{
    K <- length(len)
    unstr.row <- NULL
    unstr.col <- NULL
    ml <- max(len)
    sl2 <- sum(len^2)
    for (i in 2:ml) {
        unstr.row <- c(unstr.row, 1:(i - 1))
        unstr.col <- c(unstr.col, rep(i, each = i - 1))
    }
    unstr.row <- c(unstr.row, 1:ml)
    unstr.col <- c(unstr.col, 1:ml)
    xvec <- c(alpha.new, rep(1, ml))
    biggestMat <- forceSymmetric(sparseMatrix(i = unstr.row, 
        j = unstr.col, x = xvec))
    mat.sizes <- sort(unique(len))
    corr.vec <- vector("numeric", sl2)
    mat.inverses <- list()
    index <- vector("numeric", K + 1)
    index[1] <- 0
    index[2:K] <- (cumsum(len^2) - len^2)[2:K]
    index[K + 1] <- sl2
    for (i in 1:length(mat.sizes)) {
        tmp <- biggestMat[1:mat.sizes[i], 1:mat.sizes[i]]
        mat.inverses[[i]] <- as.vector(solve(tmp))
    }
    mat.finder <- match(len, mat.sizes)
    corr.vec <- unlist(mat.inverses[mat.finder])
    return(as(sparseMatrix(i = row.vec, j = col.vec, x = corr.vec), 
        "symmetricMatrix"))
}


getAlphaInvUser <- function (alpha.new, len, struct.vec, user.row, user.col, row.vec, 
    col.vec) 
{
    K <- length(len)
    ml <- max(len)
    sl2 <- sum(len^2)
    user.row <- c(user.row, 1:ml)
    user.col <- c(user.col, 1:ml)
    xvec <- rep.int(0, length(struct.vec))
    for (i in 1:length(alpha.new)) {
        xvec[which(struct.vec == i)] <- alpha.new[i]
    }
    xvec <- c(xvec, rep(1, ml))
    biggestMat <- forceSymmetric(sparseMatrix(i = user.row, j = user.col, 
        x = xvec))
    mat.sizes <- sort(unique(len))
    corr.vec <- vector("numeric", sl2)
    mat.inverses <- list()
    for (i in 1:length(mat.sizes)) {
        tmp <- biggestMat[1:mat.sizes[i], 1:mat.sizes[i]]
        mat.inverses[[i]] <- as.vector(solve(tmp))
    }
    mat.finder <- match(len, mat.sizes)
    corr.vec <- unlist(mat.inverses[mat.finder])
    return(as(sparseMatrix(i = row.vec, j = col.vec, x = corr.vec), 
        "symmetricMatrix"))
}


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
getSandwich <- function (YY, XX, eta, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, 
    VarFun, hessMat, StdErr, dInvLinkdEta, BlockDiag, sqrtW) 
{
    diag(dInvLinkdEta) <- InvLinkDeriv(eta)
    mu <- InvLink(eta)
    diag(StdErr) <- sqrt(1/VarFun(mu))
    scoreDiag <- Diagonal(x = YY - mu)
    BlockDiag <- scoreDiag %*% BlockDiag %*% scoreDiag
    numsand <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% 
        XX, R.alpha.inv %*% sqrtW %*% StdErr %*% BlockDiag %*% 
        StdErr %*% sqrtW %*% R.alpha.inv %*% sqrtW %*% StdErr %*% 
        dInvLinkdEta %*% XX)
    sandvar <- t(solve(hessMat, numsand))
    sandvar <- t(solve(t(hessMat), sandvar))
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


print.summary.geem2 <- function (x, ...) 
{
    Coefs <- matrix(0, nrow = length(x$coefnames), ncol = 5)
    rownames(Coefs) <- c(x$coefnames)
    colnames(Coefs) <- c("Estimates", "Model SE", "Robust SE", 
        "wald", "p")
    Coefs[, 1] <- x$beta
    Coefs[, 2] <- x$se.model
    Coefs[, 3] <- x$se.robust
    Coefs[, 4] <- x$wald.test
    Coefs[, 5] <- x$p
    print(signif(Coefs, digits = 4))
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
    cat(" Est. Scale Parameter: ", signif(x$phi, digits = 4), 
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
    naive <- is.character(object$var)
    if (!naive && any(diag(object$var) < 0)) {
        naive <- TRUE
        warning("Some elements of robust variance estimate < 0.  Reporting model based SE.")
    }
    Coefs[, 2] <- sqrt(diag(object$naiv.var))
    if (naive) {
        Coefs[, 3] <- rep(NA, length(object$beta))
    }
    else {
        Coefs[, 3] <- sqrt(diag(object$var))
    }
    if (naive) {
        Coefs[, 4] <- Coefs[, 1]/Coefs[, 2]
    }
    else {
        Coefs[, 4] <- Coefs[, 1]/Coefs[, 3]
    }
    Coefs[, 5] <- round(2 * pnorm(abs(Coefs[, 4]), lower.tail = F), 
        digits = 8)
    colnames(Coefs) <- c("Estimates", "Model SE", "Robust SE", 
        "wald", "p")
    summ <- list(beta = Coefs[, 1], se.model = Coefs[, 2], se.robust = Coefs[, 
        3], wald.test = Coefs[, 4], p = Coefs[, 5], alpha = object$alpha, 
        corr = object$corr, phi = object$phi, niter = object$niter, 
        clusz = object$clusz, coefnames = object$coefnames, weights = object$weights, 
        biggest.R.alpha = object$biggest.R.alpha)
    class(summ) <- "summary.geem2"
    return(summ)
}


updateAlphaAR <- function (YY, mu, VarFun, phi, id, len, StdErr, p, included, 
    includedlen, includedvec, allobs, sqrtW, BlockDiag, useP) 
{
    Resid <- StdErr %*% sqrtW %*% Diagonal(x = YY - mu)
    denom <- phi * (sum(band(triu(included %*% BlockDiag %*% 
        included, k = 1), k1 = 1, k2 = 1)) - useP * p)
    num <- sum(band(triu(Resid %*% BlockDiag %*% Resid, k = 1), 
        k1 = 1, k2 = 1))
    alpha <- num/denom
    return(as.numeric(alpha))
}


###Modification: updateAlphaEX is taken from geeM version 0.08 to use the variance weights
updateAlphaEX <- function (YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, 
    included, includedlen, sqrtW, useP) 
{
    Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - 
        mu)
    denom <- phi * (sum(triu(included %*% BlockDiag %*% included, 
        k = 1)) - useP * p)
    alpha <- sum(triu(Resid %*% BlockDiag %*% Resid, k = 1))
    alpha.new <- alpha/denom
    return(alpha.new)
}


updateAlphaMDEP <- function (YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, 
    m, included, includedlen, allobs, sqrtW, useP) 
{
    Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - 
        mu)
    BlockDiag <- Resid %*% BlockDiag %*% Resid
    alpha.new <- vector("numeric", m)
    for (i in 1:m) {
        if (sum(includedlen > i) > p) {
            bandmat <- drop0(band(BlockDiag, i, i))
            if (allobs) {
                alpha.new[i] <- sum(bandmat)/(phi * (sum(as.numeric(len > 
                  i) * (len - i)) - useP * p))
            }
            else {
                alpha.new[i] <- sum(bandmat)/(phi * (length(bandmat@i) - 
                  useP * p))
            }
        }
        else {
            bandmat <- drop0(band(BlockDiag, i, i))
            if (allobs) {
                alpha.new[i] <- sum(bandmat)/(phi * (sum(as.numeric(len > 
                  i) * (len - i))))
            }
            else {
                alpha.new[i] <- sum(bandmat)/(phi * length(bandmat@i))
            }
        }
    }
    return(alpha.new)
}


updateAlphaUnstruc <- function (YY, mu, VarFun, phi, id, len, StdErr, Resid, p, BlockDiag, 
    included, includedlen, allobs, sqrtW, useP) 
{
    Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - 
        mu)
    ml <- max(len)
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
    if (sum(includedlen == max(len)) <= p) {
        stop("Number of clusters of largest size is less than p.")
    }
    for (i in 1:lalph) {
        newrow <- index[which(len >= col.vec[i])] + row.vec[i]
        newcol <- index[which(len >= col.vec[i])] + col.vec[i]
        bdtmp <- BlockDiag[cbind(newrow, newcol)]
        if (allobs) {
            denom <- (phi * (length(newrow) - useP * p))
        }
        else {
            denom <- (phi * (sum(bdtmp != 0) - useP * p))
        }
        alpha.new[i] <- sum(bdtmp)/denom
    }
    return(alpha.new)
}


updateAlphaUser <- function (YY, mu, phi, id, len, StdErr, Resid, p, BlockDiag, 
    row.vec, col.vec, corr.list, included, includedlen, allobs, 
    sqrtW, useP) 
{
    Resid <- StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - 
        mu)
    ml <- max(len)
    BlockDiag <- Resid %*% BlockDiag %*% Resid
    alpha.new <- vector("numeric", length(corr.list))
    index <- cumsum(len) - len
    for (i in 1:length(alpha.new)) {
        newrow <- NULL
        newcol <- NULL
        for (j in 1:length(corr.list[[i]])) {
            newrow <- c(newrow, index[which(len >= col.vec[corr.list[[i]]][j])] + 
                row.vec[corr.list[[i]][j]])
            newcol <- c(newcol, index[which(len >= col.vec[corr.list[[i]]][j])] + 
                col.vec[corr.list[[i]][j]])
        }
        bdtmp <- BlockDiag[cbind(newrow, newcol)]
        if (allobs) {
            denom <- phi * (length(newrow) - useP * p)
        }
        else {
            denom <- phi * (sum(bdtmp != 0) - useP * p)
        }
        alpha.new[i] <- sum(bdtmp)/denom
    }
    return(alpha.new)
}


####Modification: The updateBeta function uses the variance weights like in the original geeM Version 0.08
updateBeta <- function (YY, XX, beta, off, InvLinkDeriv, InvLink, VarFun, R.alpha.inv, 
    StdErr, dInvLinkdEta, tol, sqrtW, restriction)  #restriction=list(L=NULL,r=dim(L)[1])
{
    beta.new <- beta
    conv = F
    for (i in 1:10) {
        eta <- as.vector(XX %*% beta.new) + off
        diag(dInvLinkdEta) <- InvLinkDeriv(eta)
        mu <- InvLink(eta)
        diag(StdErr) <- sqrt(1/VarFun(mu))
        hess <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% 
            XX, R.alpha.inv %*% sqrtW %*% StdErr %*% dInvLinkdEta %*% 
            XX)
        esteq <- crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% 
            XX, R.alpha.inv %*% sqrtW %*% StdErr %*% (YY - mu))
        #Modification to enable score test for the null hypotheses L%*%beta==r
        if(is.null(restriction[[1]])) {
            update <- solve(hess, esteq)
        } else {
            L<-restriction[[1]]
            r<-restriction[[2]]
            InvHess<-solve(hess)
            lambda<-solve(L%*%InvHess%*%t(L),L%*%beta.new - r + L%*%InvHess%*%esteq)
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


updatePhi <- function (YY, mu, VarFun, p, StdErr, included, includedlen, sqrtW, 
    useP) 
{
    nn <- sum(includedlen)
    resid <- diag(StdErr %*% included %*% sqrtW %*% Diagonal(x = YY - 
        mu))
    phi <- (1/(sum(included) - useP * p)) * crossprod(resid, 
        resid)
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
	W<-x$sandw.arg$dInvLinkdEta%*%x$sandw.args$StdErr %*% x$sandw.args$sqrtW %*% x$sandw.args$R.alpha.inv %*% x$sandw.args$sqrtW %*% x$sandw.args$StdErr
	if(biascorr==TRUE) {
		#Mancl and DeRouen bias corrected covariance matrix estimator
		#H Matrix fuer Mancl und DeRouen Methode (Biaskorrektur) berechnen:
		#Daten muessen dafuer nach id sortiert sein!!!
		K<-length(x$sandw.arg$len)
		Xcluster<-vector(mode="list",length=K)
		startstop<-c(0,cumsum(x$sandw.arg$len))
		for(i in 1:K) Xcluster[[i]]<-x$sandw.arg$X[(startstop[i]+1):startstop[i+1],,drop=FALSE] #wenn man eine Zeile herausholt, wird es als sonst Spaltenvektor interpretiert
		Xdiag<-bdiag(Xcluster)

		hessInv<-solve(x$sandw.args$hess)
				
		#Bias Korrektur bei restringierter Schaetzung, der slot restriction enthaelt L und r, wenn es eine Restriktion gibt, sonst ist er NULL
		#x$restricion[[1]] ist L, x$restricion[[2]] ist r
		if(!is.null(x$restricion)) {
			hessInvAdj<-hessInv-hessInv%*%t(x$restricion[[1]])%*%solve(x$restricion[[1]]%*%hessInv%*%t(x$restricion[[1]]))%*%x$restricion[[1]]%*%hessInv
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
		W<-W%*%solve(I_n - H)
	}
	W
}

getSdiag<-function(x) {
	mu <- x$sandw.arg$InvLink(x$sandw.arg$eta)
	scoreDiag <- Diagonal(x = x$sandw.arg$Y - mu)
}

getID<-function(x) x$sandw.args$id

getHess<-function(x) x$sandw.args$hess


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
#' data(lesions)
#' m1<-geem2(clearance~trt,id=id,data=lesions,family=binomial,corstr="independence")
#' m2<-geem2(pain~trt,id=id,data=lesions[lesions$lesion==1,],family=gaussian,corstr="independence")
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
	#wegen Strukrur 	0		0	1 1  0
	#			a (0, a) =	  a	1 1    a 
	Upart<-t(X.stack)%*%W.stack%*%S.stack
	B<-Upart%*%BD.stack%*%t(Upart)
	#B<-t(X.stack)%*%W.stack%*%S.stack%*%BD.stack%*%S.stack%*%t(W.stack)%*%X.stack

	#Sandwich, aeusserer Teil
	Hess.stack<-bdiag(lapply(x,getHess))
	A<-solve(Hess.stack)

	#ganzes Sandwich
	V<-A%*%B%*%A
	#Das Ergebnis ist oft im Bereich von 10^-15 asymmetrisch, sind wohl numerische Ungenauigkeiten.
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
vcov.mmmgee<-function(x) x$V
coef.mmmgee<-function(x) x$beta


getZeror<-function(L) rep(0,dim(L)[1])
getp<-function(m) length(coef(m))


mmmscore.test<-function(Modelle,L.list,r.list=NULL,tol=10^(-8),type="quadratic",alternative="undirected",biascorr=FALSE,...) {
	if(is.null(r.list)) r.list<-lapply(L.list,getZeror)
	#Fit models under restriction
	M<-length(Modelle)
	ModelleRestr<-vector(length=M,mode="list")
	for(i in 1:M) ModelleRestr[[i]]<-update(Modelle[[i]],restriction=list(L=L.list[[i]],r=r.list[[i]]),conv.criterion="difference",tol=tol)
			
	MMM<-mmmgee(ModelleRestr,biascorr=biascorr)
		
	getScore<-function(x) x$sandw.args$score[,1]
	#getScore(xRestr[[1]])
	score<-unlist(lapply(ModelleRestr,getScore))
	L<-bdiag(L.list)
	df<-dim(L)[1]
	
	#Quadratic form Score Test
	#For the robust Score Test Statistic, see e.g. D. Boos 1992, On Generalized Score Tests, The American Statistician: Vol 46, No 4
	if(type=="quadratic") {
		T<-as.numeric( t(score) %*% MMM$A%*%t(L)%*%solve( L%*%MMM$V%*%t(L) )%*%L%*%MMM$A %*% score )
		p<-1-pchisq(T,df=df) #The rank of L is checked in the geem function.
		retval<-list(test=data.frame(Chisq=T,df=df,p=p),restricted.coef=MMM$beta,score=score,type=type)	
	}
	
	#Max Score Test
	if(type=="maximum") {
		Vscore<-L%*%MMM$V%*%t(L)
		Cscore<-as.matrix(cov2cor(Vscore))
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
mmmchisq.test.simple<-function(x,L,r=NULL,approximation=c("Chisq","F","scaled.F"),df2=NA) {
	approximation<-match.arg(approximation,choices=c("Chisq","F","scaled.F"))
	T<-as.numeric( t(L%*%x$beta-r)%*%solve(L%*%x$V%*%t(L))%*%(L%*%x$beta-r) )
	df1<-dim(L)[1]
	if(approximation=="Chisq") {
		p<-1-pchisq(T,df=df1) #The rank of L is checked in the geem function.
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
		p<-1-pmvnorm(lower=lower,upper=upper,sigma=as.matrix(cov2cor(Vcontr)),...)
		retval<-list(test=data.frame(maxZ=maxz,p=p),type="MaxZ",alternative=alternative)
		names(retval$test)[1]<-c("Max|Z|","MaxZ","MinZ")[alternative==c("undirected","greater","less")]
	} else {
		p<-1-pmvt(lower=lower,upper=upper,df=floor(df2),corr=as.matrix(cov2cor(Vcontr)),...)
		retval<-list(test=data.frame(maxT=maxz,df1=df,df2=df2,p=p),type="MaxT",alternative=alternative)
		names(retval$test)[1]<-c("Max|T|","MaxT","MinT")[alternative==c("undirected","greater","less")]
	}
	retval
}





###Example Data
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
#' @usage data(lesions)
#'
#' @format A data frame.
#'
#' @keywords datasets
#'
#' @examples
#' data(lesions)
#' head(lesions)
"lesions"










######

