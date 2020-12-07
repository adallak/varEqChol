require(mgcv)
require(varband)
## Implements rowdeletion algorithm as in Davis et.al (2005)



#' Implements rowdeletion algorithm as in Davis et.al (2005)
#'
#' @param L - lower triangular Cholesky factor
#' @param k - index for row and column
#' @param const - constant value
#' @export
#' @return - Updated Cholesky factor
rowdeletion <- function(L, k, const = 0.001)
{
      p = dim(L)[1]
      new.L = L
      new.L[k,] = new.L[, k] = 0
      new.L[k, k] = const
      if(k < p)
      {
            omega = L[(k + 1):p, k]
            new.L[(k + 1):p,(k+1) : p] = t(cholup(t(L[(k + 1):p,(k+1) : p]), omega, TRUE))
      }
      return(new.L)
}


#' Returns ordering from the Cholesky factor
#' @param L - Cholesky factor of precision matrix
#' @param T - Cholesky factor of covariance matrix
#' @param p - dimensionality of the data
get_ordering_ch <- function(L, T, p)
{
      ord = c()
      for ( j in 1 : p){
            chol_i = which.min(colSums(L^2))
            ord[j] = chol_i
            T = rowdeletion(T,chol_i)
            L       = solve(T)
            ord[j] = chol_i
      }
      return(rev(ord))
}

### Returns ordering from the Cholesky factor



#' Returns ordering from the Cholesky factor
#'
#' @param X - n x p data
#' @param type - type of dimension
#' @param lamlist - A list of non-negative tuning patameters for CV
#' @param nlam    - if lamlist = NULL, create list with length equal to nlam
#' @param flmin   - If lamlist is not provided, create a lamlist with ratio of the smallest and largest
#' @param flmin  - in the list equal to flmin. Default is 0.01.
#' @param gamma  - parameter for eBIC. Active if crit = "eBIC"
#' @param refine - if TRUE, refines edges after CSCS algorithm
#' @param mtd    - parameters for refine
#' @param alpha  - parameters for refine
#' @param threshold - threshold value
#' @param FCD     - parameters for refine
#' @param precmtd - parameters for refine
#' @param crit   - criterion for penalty parameter selection
#' @export
#' @return - ordering, lower triangular Cholesky factor and adjacency matrix
eqvarDAG_ch <- function(X, type = c("low", "high"), crit = c("CV","eBIC"), gamma = 0.5, refine = FALSE, lamlist = NULL,
                        nlam = 30, flmin = 0.01, mtd="dlasso",alpha=0.05,
                        threshold=1e-1, FCD=TRUE, precmtd="sqrtlasso")
{
      p = dim(X)[2]
      n = dim(X)[1]
      ord = c()
      if (type == "low")
      {
            S = cov(X)
            T = t(chol(S))
            L = solve(T)

      }else{
            if (crit == "CV")
            {
               L = varband_cv(X, w = FALSE, lasso = TRUE, lamlist = lamlist, nlam = nlam,
                       flmin = flmin, folds = NULL, nfolds = 5)$L_refit
            }else{
               L = varband_ebic(X, w = FALSE, lasso = TRUE, lamlist = lamlist, nlam = nlam,
                                  flmin = flmin, folds = NULL, nfolds = 5, gamma = gamma)$L_refit
            }
            T = solve(L)
      }

      ord = get_ordering_ch(L, T, p)
      if(isTRUE(refine))
      {
            adj=DAG_from_Ordering(X,TO = ord,mtd,alpha,threshold,FCD,precmtd)
      }else{
         if (crit == "CV")
         {

            L = varband_cv(X[,ord], w = FALSE, lasso = TRUE, lamlist = lamlist, nlam = nlam,
                           flmin = flmin, folds = NULL, nfolds = 5)$L_refit
         }else{
            L = varband_ebic(X[, ord], w = FALSE, lasso = TRUE, lamlist = lamlist, nlam = nlam,
                             flmin = flmin, folds = NULL, nfolds = 5, gamma = gamma)$L_refit
         }

            B = diag(p) - L / (diag(L))
            adj = B * (abs(B) > threshold) != 0
      }
      ## obtain adjacency matrix
      return(list("ord" = ord, "L" = L, "adj" = adj))

}

###############
### Fit DAG using topological ordering
###############
#' Infer  DAG using topological ordering
#' @param X: data in n x p matrix
#' @param TO: topological ordering
#' @param alpha: desired selection significance level
#' @param mtd: methods for learning DAG from topological orderings.
#'  "ztest": (p<n) [Multiple Testing and Error Control in Gaussian Graphical Model Selection. Drton and Perlman.2007]
#'  "rls": (p<n) fit recursive least squares using ggm package and threshold the regression coefs
#'  "chol": (p<n) perform cholesky decomposition and threshold the regression coefs
#'  "dlasso": debiased lasso (default with FCD=True and precmtd="sqrtlasso");
#'   "lasso": lasso with fixed lambda from [Penalized likelihood methods for estimation of sparse high-dimensional directed acyclic graphs. Shojaie and Michailidis. 2010];
#'   "adalasso": adaptive lasso with fixed lambda from [Shojaie and Michailidis. 2010];
#'   "cvlasso": cross-validated lasso from glmnet;
#'    "scallasso": scaled lasso.
#' @param threshold: only used in rls and chol. the hard threshold level.
#' @param FCD: only used in debiased lasso,  the FCD procedure [False Discovery Rate Control via Debiased Lasso. Javanmard and Montanari. 2018]
#' or use individual tests to select support.
#' @param precmtd: only used in debiased lasso, how to compute debiasing matrix
#'               "cv": node-wise lasso w/ joint 10 fold cv
#'               "sqrtlasso": square-root lasso(no tune, default)
#' @return Adjacency matrix with ADJ[i,j]!=0 iff i->j
DAG_from_Ordering<-function(X,TO,mtd="ztest",alpha=0.05,
                            threshold=1e-1,FCD=NULL,precmtd=NULL){
      n=dim(X)[1]
      p=dim(X)[2]
      if (p!=length(TO)){stop("length mismatch")}
      if (mtd=="ztest"){
            # sidak
            C=cor(X)
            adj=matrix(0,p,p)
            for (i in 2:p){
                  u=TO[i]
                  for (j in 1:(i-1)){
                        v = TO[j]
                        s = setdiff(TO[seq(i-1)],v)
                        pval = 1-(2*pnorm(abs(pcalg::zStat(u,v,s,C=C,n=n)))-1)^(p*(p-1)/2)
                        adj[v,u]=ifelse(pval<alpha,1,0)
                  }
            }
            return(adj!=0)
      }
      if (mtd=="chol"){
            Sigma=cov(X)
            B = solve(chol(Sigma[TO,TO])[order(TO),order(TO)])
            gm = diag(p)-B%*%diag(1/diag(B))
            return(gm*(abs(gm)>threshold)!=0)
      }
      if (mtd=="rls"){
            gm = upper.tri(matrix(0,p,p))[order(TO),order(TO)]
            colnames(gm)=rownames(gm)=colnames(X)
            return(abs(t(ggm::fitDag(gm,cov(X),dim(X)[1])$Ahat))-diag(p)>threshold)
      } else {
            # dblasso
            if (is.null(FCD)){FCD="T"}
            if (is.null(precmtd)){precmtd="sqrtlasso"}
            gm = matrix(0,p,p)
            gm[TO[1],TO[2]]=anova(lm(X[,TO[2]]~X[,TO[1]]))$`Pr(>F)`[1]<alpha
            if(p==2){return(gm)}
            for (i in 3:p){
                  gm[TO[1:(i-1)],TO[i]]=
                        vselect(X[,TO[1:i-1]],X[,TO[i]],alpha=alpha,p_total = p,
                                selmtd = mtd,FCD = FCD,precmtd = precmtd)$selected
            }
            return(gm!=0)
      }
}

######################################
#' Compares the estimated and true adjacency matrices
#' @param estAdj - estimated adjacency matrix
#' @param trueAdj - true adjacency matrix
#'
#' @return True Positive, False Positive and True Discovery rates
#'
#' @export
compareGraph <- function (estAdj, trueAdj)
{

      #### This function compares esitmated
      ### adjacency matrix of A with the true Adj matrix
      ### of coefficient matrix A
      ml <- estAdj
      mt <- trueAdj
      p <- dim(ml)[2]
      mt[mt != 0] <- rep(1, sum(mt != 0))
      ml[ml != 0] <- rep(1, sum(ml != 0))
      diffm <- ml - mt
      nmbTrueGaps <- (sum(mt == 0) - p)/2
      fpr <- if (nmbTrueGaps == 0)
            1
      else (sum(diffm > 0)/2)/nmbTrueGaps
      diffm2 <- mt - ml
      nmbTrueEdges <- (sum(mt == 1)/2)
      tpr <- if (nmbTrueEdges == 0)
            0
      else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
      trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
      tdr <- if (sum(ml == 1) == 0) {
            if (trueEstEdges == 0)
                  1
            else 0
      }
      else trueEstEdges/(sum(ml == 1)/2)
      return(list(tpr = tpr, fpr = fpr, tdr = tdr))
}


##########################################################
varband_ebic <- function(x, w = FALSE, lasso = TRUE, lamlist = NULL, nlam = 40,
                 flmin = 1e-2, folds = NULL, nfolds = 5, gamma = 0.5)
{
   n <- nrow(x)
   p <- ncol(x)

   S <- crossprod(scale(x, center=TRUE, scale=FALSE)) / n
   if(is.null(folds))
      folds <- makefolds(n, nfolds = nfolds)
   nfolds <- length(folds)

   if (is.null(lamlist)) {
      lam_max <- lammax(S = S)
      lamlist <- pathGen(nlam = nlam, lam_max = lam_max,
                         flmin = flmin, S = S)
   } else {
      nlam <- length(lamlist)
   }

   errs_fit <- matrix(NA, nlam, nfolds)
   errs_refit <- matrix(NA, nlam, nfolds)

   # error function is the negative log Gaussian likelihood
   for (i in seq(nfolds)) {
      # train on all but i-th fold:
      x_tr <- x[-folds[[i]],]
      meanx <- colMeans(x_tr)
      x_tr <- scale(x_tr, center = meanx, scale = FALSE)
      S_tr <- crossprod(x_tr) / (dim(x_tr)[1])

      path_fit <- varband_path(S = S_tr, w = w, lasso = lasso,
                               lamlist = lamlist)$path
      path_refit <- refit_path(S = S_tr, path = path_fit)

      # evaluate this on left-out fold:
      x_te <- x[folds[[i]], ]
      x_te <- scale(x_te, center = meanx, scale = FALSE)
      S_te <- crossprod(x_te) / (dim(x_te)[1])

      for (j in seq(nlam)) {
         errs_fit[j, i] <- eBIC(path_fit[, , j], S_te, n, gamma = gamma)
         errs_refit[j, i] <- eBIC(path_fit[, , j], S_te, n, gamma = gamma)
      }
   }

   m_fit <- rowMeans(errs_fit)
   se_fit <- apply(errs_fit, 1, sd) / sqrt(nfolds)
   m_refit <- rowMeans(errs_refit)
   se_refit <- apply(errs_refit, 1, sd) / sqrt(nfolds)
   ibest_fit <- which.min(m_fit)
   i1se_fit <- min(which(m_fit < m_fit[ibest_fit] + se_fit[ibest_fit]))
   ibest_refit <- which.min(m_refit)
   i1se_refit <- min(which(m_refit < m_refit[ibest_refit] + se_refit[ibest_refit]))


   fit_ebic <- varband(S = S, lambda = lamlist[ibest_fit], init = path_fit[, , ibest_fit], w = w, lasso = lasso)
   refit_ebic <- varband(S = S, lambda = lamlist[ibest_refit], init = path_refit[, , ibest_refit], w = w, lasso = lasso)
   refit_ebic <- refit_matrix(S = S, mat = refit_ebic)

   return(list(errs_fit = errs_fit, errs_refit = errs_refit,
               folds = folds, lamlist = lamlist,
               ibest_fit = ibest_fit, ibest_refit = ibest_refit,
               i1se_fit = i1se_fit, i1se_refit = i1se_refit,
               L_fit = fit_ebic, L_refit = refit_ebic))
}

makefolds <- function(n, nfolds) {
   # divides the indices 1:n into nfolds random folds of about the same size.
   nn <- round(n / nfolds)
   sizes <- rep(nn, nfolds)
   sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
   b <- c(0, cumsum(sizes))
   ii <- sample(n)
   folds <- list()
   for (i in seq(nfolds))
      folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
   folds
}

eBIC <- function (L, S, n, lambda, gamma = 0.5){
   # Calculate the negative log-Gaussian likelihood with
   # precision matrix Omega and sample covariance S
   p = dim(L)[1]
   s = sum(which(abs(L) > 0)) - p
   return(-2 * sum(log(diag(L))) + sum(diag(crossprod(L) %*% S))
          + s * log(n) + 4 * s * gamma * log(p))
}

lammax <- function(S){
   #### This function calculates the max value in the tuning parameter list
   # such that the estimator L_{\lambda} is a diagonal matrix
   # NOTE: this is not necessarily true, but generally
   # a upper bound of the value we are looking for.

   # Args:
   #     S: the p-by-p sample covariance matrix

   p <- ncol(S)
   sighat <- rep(NA, p-1)
   for (r in seq(2, p)){
      sighat[r-1] <- max(abs(S[(1:(r-1)), r]))/sqrt(S[r, r])
   }
   2 * max(sighat)
}

pathGen <- function(nlam, lam_max, flmin, S){
   # Generate a path of lambda, with
   # nlam/2 decreasing exponentially
   # nlam/2 decreasing linearly
   # lam_max <- lammax(S)
   lamlist_lin <- lam_max * exp(seq(0, log(flmin), length = nlam/2))
   lamlist_exp <- seq(lam_max - 1e-8, lam_max*flmin - 1e-8, length.out = nlam/2)
   return(sort(unique(c(lamlist_lin, lamlist_exp)), decreasing = TRUE))
}

THRESH <- 1e-10
# #' Refit a row problem with given support
# #'
# #' @param r row index
# #' @param ind the index set of support
# #' @param S p-by-p sample covariance matrix
# #' @param delta a nonnegative tuning parameter for the ridge penalty for numerical stability, default is 0.01
refit_row <- function(r, ind, S, delta = 0.01){
   # Refitting by plain MLE with
   # bandwidth k = r-1-J given our estimator
   p <- ncol(S)
   res <- rep(0, r)
   if(is.null(ind)){
      res[r] <- 1/sqrt(S[r, r])
   }
   else{
      # if J < r-1
      # ind <- (J+1):(r-1)
      # If S[ind,ind] is not invertible
      # add a little ridge penalty to that
      tmpVec <- solve((S[ind, ind] +
                          delta * diag(rep(1, length(ind)))),
                      S[ind, r])
      res[r] <- 1/sqrt(S[r, r] + delta -
                          crossprod(S[r, ind], tmpVec))
      res[ind] <- -tmpVec * res[r]
   }
   return(res)
}

# #' Refit the estimate of lower triangular matrix L with given support
# #'
# #' @param S p-by-p sample covariance matrix
# #' @param mat p-by-p estimate of lower triangular matrix L
# #' @param delta a nonnegative tuning parameter for the ridge penalty for numerical stability, default is 0.01
refit_matrix <- function(S, mat, delta = 0.01){
   p <- ncol(S)
   refit <- matrix(0, p, p)
   refit[1, 1] <- 1/sqrt(S[1, 1])
   for(r in seq(2, p)){
      ind <- c()
      for(j in seq(1,r-1)){
         if(abs(mat[r, j]) >= THRESH){
            ind <- c(ind, j)
         }
      }
      refit[r, 1:r] <- refit_row(r = r, ind = ind,
                                 S = S, delta = delta)
   }
   return(refit)
}

# #' Refit a path of estimates of lower triangular matrix L with given support
# #'
# #' @param  S: p-by-p sample covariance matrix
# #' @param  path: a list of p-by-p estimate of lower triangular matrix L along a path of tuning parameters
# #' @param delta: a nonnegative tuning parameter for the ridge penalty for numerical stability, default is 0.01
refit_path <- function(S, path, delta = 0.01){
   p <- dim(path)[1]
   nlam <- dim(path)[3]
   refit <- array(NA, c(p, p, nlam))
   for(i in seq(nlam)){
      refit[, , i] <- refit_matrix(S, path[, , i], delta)
   }
   return(refit)
}


hammingDistance <- function(G1,G2)
   # hammingDistance(G1,G2)
   #
   # Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1!!!!
   #
   # INPUT:  G1, G2     adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
   #
   # OUTPUT: hammingDis Hamming Distance between G1 and G2
   #
   # Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]
   # All rights reserved.  See the file COPYING for license terms.
{
   allMistakesOne <- FALSE
   if(allMistakesOne)
   {
      Gtmp <- (G1+G2)%%2
      Gtmp <- Gtmp + t(Gtmp)
      nrReversals <- sum(Gtmp == 2)/2
      nrInclDel <- sum(Gtmp == 1)/2
      hammingDis <- nrReversals + nrInclDel
   } else
   {
      hammingDis <- sum(abs(G1 - G2))
      # correction: dist(-,.) = 1, not 2
      hammingDis <- hammingDis - 0.5*sum(G1 * t(G1) * (1-G2) * t(1-G2) + G2 * t(G2) * (1-G1) * t(1-G1))
   }
   return(hammingDis)
}

# ER graph
randomDAG2_er <- function(p,probConnect)
{
   # This function is modified from randomDAG2 function by Jonas Peters
   DAG <- diag(rep(0,p))
   causalOrder <- sample(p)
   for(i in 3:(p))
   {
      node <- causalOrder[i]
      possibleParents <- causalOrder[1:(i-1)]
      Parents <- possibleParents[rbinom(length(possibleParents),1,probConnect)==1]
      DAG[Parents,node] <- rep(1,length(Parents))
   }
   node <- causalOrder[p-1]
   ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
   DAG[causalOrder[1],causalOrder[2]] <- 1

   return(list(DAG=DAG,TO=causalOrder))
}
# Chain graph
randomDAG2_chain <- function(p,probConnect)
{
   # This function is modified from randomDAG2 function by Jonas Peters
   DAG <- diag(rep(0,p))
   causalOrder <- sample(p)
   DAG[causalOrder[1],causalOrder[2]] <- 1
   for(i in 3:(p))
   {
      node <- causalOrder[i]
      possibleParents <- causalOrder[1:(i-1)]
      possibleParents <- possibleParents[which(rowSums(DAG[possibleParents,])<4)]
      if (length(possibleParents)>0){
         Parents <- sample(possibleParents,min(length(possibleParents),2))
         DAG[Parents,node] <- rep(1,length(Parents))
      }
      DAG[causalOrder[i-1],node]<-1
   }
   return(list(DAG=DAG,TO=causalOrder))
}
# Hub-and-chain graph
randomDAG2_hub <- function(p,probConnect)
{
   # This function is modified from randomDAG2 function by Jonas Peters
   DAG <- diag(rep(0,p))
   causalOrder <- sample(p)
   Z<-10
   for(i in 1:(p))
   {
      node <- causalOrder[i]
      DAG[causalOrder[i-1],node]<-1
      if (i>2){DAG[causalOrder[sample(seq(min(i-1,Z)),2)],node]<-1}
   }
   DAG[causalOrder[1],causalOrder[2]] <- 1
   return(list(DAG=DAG,TO=causalOrder))
}

###############
### Generate data
###############
Bmin<-0.5
get_DAGdata<-function(n,p,pc,type='hub',err='g'){
   if (type=='hub'){
      D<-randomDAG2_hub(p,pc)
   } else if (type=='chain') {
      D<-randomDAG2_chain(p,pc)
   } else {
      D<-randomDAG2_er(p,pc)
   }
   truth<-D$DAG
   TO<-D$TO
   if (err == 'nong'){
      errs <- matrix((rbinom(p * n,1,0.5)*2-1)*sqrt(0.8), nrow = p, ncol = n)
   } else {
      errs <- matrix(rnorm(p * n), nrow = p, ncol = n)
   }
   B<-t(truth)
   B[B==1]<-runif(sum(truth),Bmin,1)*(2*rbinom(sum(truth),1,0.5)-1)
   X <- solve(diag(rep(1, p)) - B, errs)
   X <- t(X)
   return(list(truth=truth,B=B,X=X,TO=TO))
}
