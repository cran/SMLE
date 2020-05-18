#'Joint feature screening via sparse maximum likelihood estimation for GLMs
#'
#' Input a \eqn{n \times 1} response Y and a \eqn{n \times p} feature matrix; the function uses SMLE to
#' retain only a set of \eqn{k<n} features that seem to be most relevant for a GLM. It thus serves as a pre
#' processing step for an elaborative analysis.In SMLE, the joint effects between features are naturally
#' accounted; this makes the screening more reliable. The function uses the efficient iterative hard
#' thresholding (IHT) algorithm with step parameter adaptively tuned for fast convergence.
#' Users can choose to further conduct an elaborative selection after SMLE-screening. See
#' smle_select for more details.
#'
#' @import glmnet mnormt
#' @importFrom glmnet glmnet
#' @importFrom grDevices dev.new dev.off rainbow
#' @importFrom graphics barplot legend lines par plot title
#' @importFrom stats dpois glm lm model.frame model.matrix quantile rbinom rnorm rpois
#' @importFrom mnormt rmnorm
#'
#' @details
#' When given the input data X and Y, SMLE first runs glmnet to find a initial values to start the Iterative
#' hard thresholding algrithm(IHT). The step size in IHT is selected by an adaptively algorithm named Ucheck. For more information about IHT please see the reference given below.
#'
#' In Ucheck, the algorithm initializes a starting step size \code{U} by single value decomposition of Retained data X
#' or approximate it by rowsum. Then, in each iteration, we search the largest step by compromising the step size at a given rate to
#' ensure that likelihood is increasing at the direction.
#'
#' The algorithm terminates when the difference between two successive estimated coefficient vectors are within the tolerance. Because the dimension of coefficient vectors depends
#' on the number of retained features \code{k}, we can accelerate the convergence by using a heuristic stopping criteria, that is \code{fast=TRUE }. The trade off of early stopping is
#' that convergence can not be theoretically guaranteed.
#'
#' For categorical features, this function first converts categorical columns into several numerical columns by dummy coding(\code{codeingtype='standard'}) or devation
#' coding (\code{codingtype='DV'}), and then, measures importance based on the effect size of those numercial columns. Users can treat them
#' as a group(\code{group = TRUE}) or individually by \code{penalize_mod}. When \code{group = TRUE}, Categorical features' coefficients are
#' computed by averaging all the dummy variable. \deqn{\beta_i = \frac{1}{\sqrt{J}}*\sqrt{\beta_1^2+..+\beta_J^2}} where J is the level of i-th categorical data.
#'
#' Users may also wish to further select important features after screening. This can be done either by using \code{selection=TRUE} or by applying the smle_select function on the output of SMLE.
#'
#' @param Y The response vector of dimension \eqn{n \times 1}. Quantitative for
#' \code{family ='gaussian'}, non-negative counts for \code{family ='poisson'},
#' binary (0-1) for \code{family ='binomial'}. Input Y should be \code{'numeric'}.
#' @param X The design matrix with dimensions n and p. Each row is an observation vector. Each column is a covariate. The input should be
#' the object of "matrix" for numerical data, and "data.frame" for categorical data (or a mixture of numerical andcategorical data).
#' The algorithm will treat covariates having class "factor" as categorical data and extend the data frame dimension by many columns as needed for the coding type and levels of the factor.
#' @param k Number of features to be retained after screening. Default is \eqn{\frac{1}{2}\log(n)n^{1/3}}
#' @param family Model assumption between Y and X; the default model is Gaussian linear.
#' @param categorical Logical flag whether the input design matrix includes categorical covariates, if \code{categorical= TRUE},
#' intercept is added to the model. Default is NULL.
#' @param U_rate Decreasing rate in Ucheck. More details see Ucheck
#' @param keyset A vector to indicate a set of key features that do not participate in feature screening and are forced to remain in the model. Default is null.
#' @param intercept A vector to indicate whether toan intercept be used in the model. An intercept will not participate in screening.
#' @param group Logical flag for whether to treat categorical variable as a group. (Only for categorical data, see details). Default is TRUE.
#' @param codingtype Coding types for categorical variable, Default is "DV".
#' \code{Codingtype = "all"} Convert each level to a 0-1 vector.
#' \code{Codingtype = "DV"} Compares each level to the grand mean.
#' \code{Codingtype = "standard"} Compares each level to the reference level.(Only for categorical data)
#' @param penalize_mod A logical flag to indicate whether adjustment is used in ranking groups of features. This augment is applicable only when cate=T with group=T; the default is true: a factor of sqrt(J) is divided from the L2 effect of a group with J members.
#' @param standardize Logical flag for feature standardization, prior to
#' performing (iterative) feature screening.  The resulting coefficients are
#' always returned on the original scale. Default is \code{standardize=TRUE}.
#' If variables are in the same units already, you might not wish to
#' standardize.
#' @param fast Set to TRUE to speed up the screening Default is FALSE, see details.
#' @param maxit Maximum number of iteration steps. Default is 500. Setmaxit= NULL to loosen this stopping criteria of steps.
#' @param tol Tolerance of squared sum of differences between coefficient estimates for successive iterations, when value is below tol, iterations are terminated.
#' Default is 10^(-3).Setting maxit= NULL to loosen the stopping criteria of difference between two steps.
#' @param selection A logical flag to indicate whether an elaborate selection is to be conducted by smle_select after screening.Default is FALSE.
#' @param ... Other parameter for smle_selection.
#'
#' @references
#' Xu, C. and Chen, J. (2014). The Sparse MLE for Ultrahigh-Dimensional Feature Screening, \emph{Journal of the American Statistical Association}
#' @return
#'
#' Returns a \code{"smle"} object with
#' \item{I}{Preprocessed input.
#' \code{'Y'}: Same as input Y.
#' \code{'CM'}: Design matrix of class 'matrix' for numeric (or data.frame with categorical features).
#' \code{'DM'}: A matrix with dummy variable columns added.(only if there are categorical features).
#' \code{'IM'}: Iteration Matrix for algorithm.
#' \code{'nlevel'}: number of levels for all categorical factors.
#' \code{'CI'}: indices of categorical factors in 'CM'.
#' \code{'Beta0'}: Inital value for regreesion coefficients.
#' \code{'DFI'}: indices of categorical factors in 'IM'.
#' \code{'codingtype'}: Same as input.
#'  }
#' \item{Retained_Feature_IDs}{A vector indicating the features retained after SMLE screening. The output includes both features retained by SMLE and the features forced into the final model by the users}
#' \item{Coefficients_of_Retained_Features}{The vector of coefficients for the selected features}
#' \item{Retained_Features_path}{A matrix of dimension P by Steps. Each column gives a regression coefficient over the iterations.}
#' \item{Number_of_Retained_Features}{Number of featrues after screening.}
#' \item{Intercept}{The value, if intercept = TRUE.}
#' \item{steps}{Number iterations.}
#' \item{LH}{List giving log likelihood at each step.}
#' \item{Uchecks}{Number of 'U_check's at each step.}
#'
#'
#'
#'
#'
#' @export
#'
#' @examples
#'
#' #Example
#' set.seed(123.456)
#' Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
#' Data
#' fit<-SMLE(Data$Y, Data$X, k=9, family = "gaussian")
#' fit
#' ## The important features we missed:
#' setdiff(Data$index,fit$Retained_Feature_IDs)
#' ## Check if the important featrues are retained.
#' Data$index %in% fit$Retained_Feature_IDs
#' plot(fit)
#'
#'
SMLE<-function(Y,X,k=NULL,family=c("gaussian","binomial","poisson"),
               categorical=NULL,keyset = NULL, intercept = TRUE,group= TRUE,codingtype=NULL,
               maxit = 500, tol = 10^(-4), selection=F, standardize =FALSE,fast=FALSE, U_rate=0.5,
               penalize_mod=TRUE,...){
  #-------Input preprocess-------
  family<-match.arg(family)
  call<-match.call()
  if(is.null(k)){
    k<-floor(1/2*log(dim(X)[1])*dim(X)[1]^(1/3))
  }
  if(is.null(categorical)){categorical= FALSE}
  if(standardize== TRUE){X<-Standardize(X)}
  #------------------------------
  if(sum(is.na(c(Y,X)))==TRUE){stop("NA in X or Y")}

  #-------Run Algoriathm------
  if(categorical == TRUE)
  {
    if(is.null(codingtype)){codingtype<-"DV"}
    if(!codingtype%in%c("all","standard","DV")){stop("Codingtype only in all,standard,DV")}
    fit<- ctg_fit(Y,X,k,family,categorical,keyset,maxit,tol,intercept,group,codingtype,penalize_mod)
    return(fit)
  }else{
    fit<-GPBnet(Y,X,k,family,keyset,intercept,maxit,tol,fast,U_rate)
  }
  #------Adding method-----
  fit$call=call
  class(fit)="smle"
  #-----Stage II ---------
  if(selection==F){
    return(fit)
    }else{
  return(selection(fit,...))
  }
}
