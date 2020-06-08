#'Joint feature screening via sparse maximum likelihood estimation for GLMs
#'
#'
#' Input a \eqn{n \times 1} response Y and a \eqn{n \times p} feature matrix X;
#' the function uses SMLE to retain only a set of \eqn{k<n} features that seem
#' to be most relevant for a GLM. It thus serves as a pre-processing step for an
#' elaborative analysis. In SMLE, the joint effects between features are naturally
#' accounted; this makes the screening more reliable. The function uses the
#' efficient iterative hard thresholding (IHT) algorithm with step parameter
#' adaptively tuned for fast convergence. Users can choose to further conduct
#' an elaborative selection after SMLE-screening. See \code{smle_select} for more details.
#'
#' @import glmnet
#' @importFrom glmnet glmnet
#' @importFrom grDevices dev.new dev.off rainbow
#' @importFrom graphics barplot legend lines par plot title
#' @importFrom stats dpois glm lm model.frame model.matrix quantile rbinom rnorm rpois
#' @details
#'
#' With the input Y and X, \code{SMLE} conducts joint feature screening by running
#' iterative hard thresholding algorithm (IHT), where the initial value is set to
#' be the Lasso estimate with the sparsity closest to the sample size minus one.
#'
#' In \code{SMLE}, the step parameter \eqn{u^{-1}} in IHT is adaptively tuned in
#' the same way as described in Xu and Chen (2014). Specifically, at each step,
#' we set the initial \code{u} as the max row sum of X and recursively decrease
#' the value of \eqn{u^{-1}} by \code{U_rate} to guarantee the likelihood increment.
#'
#' \code{SMLE} terminates IHT iterations when either \code{tol} or \code{maxit} is
#' satisfied. When \code{fast=TRUE}, the algorithm also stops when the non-zero
#' members of the coefficient estimates remain the same for \eqn{1_0} successive
#' iterations.
#'
#' In \code{SMLE}, categorical features are coded by dummy covariates with the
#' method specified in \code{codingtype}. Users can use \code{group} to specify
#' whether to treat those dummy covariates as a single group feature or as
#' individual features.
#' When \code{group=True} with \code{penalize_mod=True}, the effect for a group
#' of \eqn{J} dummy covariates is computed by
#'
#' \deqn{ \beta_i = \frac{1}{\sqrt{J}} \cdot \sqrt{(\beta_1)^2+...+(\beta_J)^2}}
#'
#' which will be treated as a single feature in IHT iterations.
#'
#' Since feature screening is usually a preprocessing step, users may wish to
#' further conduct an elaborative feature selection after screening. This can
#' be done by setting \code{selection=True} in SMLE or applying any existing
#' selection method on the output of \code{SMLE}.
#'
#'
#'
#' @param Y The response vector of dimension \eqn{n \times 1}. Quantitative for
#' \code{family ='gaussian'}, non-negative counts for \code{family ='poisson'},
#' binary (0-1) for \code{family ='binomial'}. Input Y should be \code{'numeric'}.
#'
#' @param X The \eqn{n \times p} feature matrix with each column denoting a feature
#' (covariate) and each row denoting an observation vector. The input should be
#' the object of "matrix" for numerical data, and "data.frame" for categorical
#' data (or a mixture of numerical and categorical data). The algorithm will
#' treat covariates having class "factor" as categorical data and extend the data
#' frame dimension by the dummy columns needed for coding the categorical features.
#'
#' @param k Number of features to be retained after screening. Default is
#' \eqn{\frac{1}{2}\log(n)n^{1/3}}
#'
#' @param family Model assumption between Y and X; the default model is Gaussian
#' linear.
#'
#' @param categorical Logical flag whether the input feature matrix includes
#' categorical features. If \code{categorical= TRUE}, a model intercept will
#' be used in the screening process. Default is NULL.
#'
#' @param U_rate Decreasing rate in tuning step parameter \eqn{u^{-1}} in IHT
#' algorithm. See details.
#'
#' @param keyset A vector to indicate a set of key features that do not
#' participate in feature screening and are forced to remain in the model.
#' Default is null.
#'
#' @param intercept A vector to indicate whether to an intercept be used in
#' the model. An intercept will not participate in screening.
#'
#' @param group Logical flag for whether to treat the dummy covariates of a
#' categorical feature as a group. (Only for categorical data, see details).
#' Default is TRUE.
#'
#' @param codingtype Coding types for categorical features; default is "DV".
#' \code{Codingtype = "all"} Convert each level to a 0-1 vector.
#' \code{Codingtype = "DV"} conducts deviation coding for each level in
#' comparison with the grant mean.
#' \code{Codingtype = "standard"} conducts standard dummy coding for each level
#' in comparison with the reference level (first level).
#'
#' @param penalize_mod A logical flag to indicate whether adjustment is used in
#' ranking groups of features. This augment is applicable only when
#' \code{categorical= TRUE} with \code{group=T}; the default is true:
#' a factor of \eqn{sqrt{J}} is divided from the \eqn{L_2} effect of a group with J members.
#'
#' @param standardize Logical flag for feature standardization, prior to
#' performing (iterative) feature screening.  The resulting coefficients are
#' always returned on the original scale. Default is \code{standardize=TRUE}.
#' If features are in the same units already, you might not wish to
#' standardize.
#'
#' @param fast Set to TRUE to enable early stop for SMLE-screening. It may help
#' to boost the screening efficiency with a little sacrifice of accuracy.
#' Default is FALSE, see details.
#'
#' @param maxit Maximum number of iteration steps. Default is 500. Set
#' \code{maxit= NULL} to loosen this protective stopping criterion.
#'
#' @param tol A tolerance level to stop the iteration, when the squared sum of
#' differences between two successive coefficient updates is below it.
#' Default is \eqn{10^{-3}}. Set \code{tol= NULL} to loosen this stopping criterion.
#'
#' @param selection A logical flag to indicate whether an elaborate selection
#' is to be conducted by \code{smle_select} after screening. Default is FALSE.
#'
#' @param ... Other parameters for \code{smle_select}
#'
#' @references
#' Xu, C. and Chen, J. (2014). The Sparse MLE for Ultrahigh-Dimensional Feature
#' Screening, \emph{Journal of the American Statistical Association}
#' UCLA: Statistical Consulting Group. \emph{CODING SYSTEMS FOR CATEGORICAL VARIABLES IN REGRESSION ANALYSIS}
#' \url{https://stats.idre.ucla.edu/spss/faq/coding-systems-for-categorical-variables-in-regression-analysis-2/}
#'
#'
#'
#' @return
#'
#'
#'
#' Returns a '\code{smle}' object with
#' \item{I}{A list of iteration information.
#'
#' \code{Y}: Same as input Y.
#'
#' \code{CM}: Design matrix of class \code{matrix} for numeric features (or data.frame with categorical features).
#'
#' \code{DM}: A matrix with dummy variable featrues added. (only if there are categorical features).
#'
#' \code{IM}: Iteration path matrix with columns recording IHT coefficient updates.
#'
#' \code{nlevel}: Number of levels for all categorical features.
#'
#' \code{CI}: Indices of categorical features in \code{CM}.
#'
#' \code{Beta0}: Inital value of regression coefficient for IHT.
#'
#' \code{DFI}: Indices of categorical features in \code{IM}.
#'
#' \code{codingtype}: Same as input.
#'  }
#'
#' \item{ID_Retained}{A vector indicating the features retained after SMLE screening.
#' The output includes both features retained by SMLE and the features specified in \code{keyset}.}
#'
#' \item{Coef_Retained}{The vector of coefficients for the retained features.}
#'
#' \item{Path_Retained}{Iteration path matrix with columns recording the coefficient updates over the IHT procedure.}
#'
#' \item{Num_Retained}{Number of retained featrues after screening.}
#'
#' \item{Intercept}{The value, if Intercept = TRUE.}
#'
#' \item{steps}{Number of iterations.}
#'
#' \item{LH}{A list of log-likelihood updates over the IHT iterations }
#'
#' \item{Uchecks}{Number of times in searching a proper \eqn{u^{-1}} at each step over the IHT iterations.}
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
#' setdiff(Data$index,fit$ID_Retained)
#' ## Check if the important featrues are retained.
#' Data$index %in% fit$ID_Retained
#' plot(fit)
#'
#'
SMLE<-function(Y , X , k=NULL , family=c("gaussian","binomial","poisson"),
               categorical = NULL , keyset = NULL, intercept = TRUE ,
               group = TRUE , codingtype = NULL ,
               maxit = 50 , tol = 10^(-4) , selection = F ,
               standardize = FALSE , fast = FALSE , U_rate=0.5 ,
               penalize_mod = TRUE , ...){

  #-------Input preprocess-------

  family<-match.arg(family)

  call<-match.call()

  if(is.null(k)){
    k<-floor(1/2*log(dim(X)[1])*dim(X)[1]^(1/3))
  }

  if(is.null(categorical)){
    Ci<-(1:dim(X)[2])[sapply(X,is.factor)]
    if( any(Ci) ){
      categorical = TRUE
    }else{
        categorical= FALSE
        }
    }

  if(standardize== TRUE){X<-Standardize(X)}
  #------------------------------
  if(sum(is.na(c(Y,X)))==TRUE){stop("NA in X or Y")}

  #-------Run Algoriathm------




  if(categorical == TRUE)

  {

    if(is.null( codingtype )){codingtype<-"DV"}

    if(!codingtype%in% c("all","standard","DV")){stop("Codingtype only in all,standard,DV")}

    fit<- ctg_fit(Y,X,k,family,categorical,keyset,maxit,tol,intercept,group,codingtype,penalize_mod,call)

    return(fit)

  }
  else{

    fit<-GPBnet(Y,X,k,family,keyset,intercept,maxit,tol,fast,U_rate)

  }

  #------Adding method-----

  fit$call=call

  class(fit)="smle"

  #-----Stage II ---------
  if(selection==F){

    return(fit)

    }else{

  return(selection(fit,...)

         )
  }
}
