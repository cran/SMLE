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
#' an elaborative selection after SMLE-screening. See \code{smle_select()} for more details.
#'
#' @importFrom glmnet glmnet
#' @importFrom grDevices dev.new dev.off rainbow
#' @importFrom graphics barplot legend lines par plot title
#' @importFrom stats dpois glm lm model.frame model.matrix quantile rbinom rnorm rpois sd .getXlevels model.response
#' @importFrom utils tail
#' @details
#'
#' With the input Y and X, \code{SMLE} conducts joint feature screening by running
#' iterative hard thresholding algorithm (IHT), where the default initial value is set to
#' be the Lasso estimate with the sparsity closest to the sample size minus one.
#'
#' In \code{SMLE}, the initial value for parameter \eqn{u^{-1}} is set to be
#' \eqn{\frac{1}{||X||^2_{\infty}}} and recursively decrease the value of 
#' \eqn{u^{-1}} by \code{U_rate} to guarantee the likelihood increment.
#'
#' \code{SMLE} terminates IHT iterations when either \code{tol} or \code{max_iter} is
#' satisfied. When \code{fast=TRUE}, the algorithm also stops when the non-zero
#' members of the coefficient estimates remain the same for 10 successive
#' iterations or the log-likelihood difference between coefficient estimates is less
#' than 0.01 * the log-likelihood increase of the first step, or the 
#' \eqn{\sqrt{k}\cdot}\code{tol} is satisfied.
#'
#' In \code{SMLE}, categorical features are coded by dummy covariates with the
#' method specified in \code{codingtype}. Users can use \code{group} to specify
#' whether to treat those dummy covariates as a single group feature or as
#' individual features.
#' When \code{group=TRUE} with \code{penalize_mod=TRUE}, the effect for a group
#' of \eqn{J} dummy covariates is computed by
#'
#' \deqn{ \beta_i = \frac{1}{\sqrt{J}} \cdot \sqrt{(\beta_1)^2+...+(\beta_J)^2}}
#'
#' which will be treated as a single feature in IHT iterations.
#'
#' Since feature screening is usually a preprocessing step, users may wish to
#' further conduct an elaborative feature selection after screening. This can
#' be done by setting \code{selection=TRUE} in \code{SMLE()} or applying any existing
#' selection method on the output of \code{SMLE()}.
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
#' @param k Total number of features (including \code{'keyset'}) to be retained 
#' after screening. Default is the largest integer not 
#' exceeding  \eqn{\frac{1}{2}\log(n)n^{1/3}}.
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
#' @param codingtype Coding types for categorical features; default is \code{'DV'}.
#' \code{Codingtype = "all"} Convert each level to a 0-1 vector.
#' \code{Codingtype = "DV"} conducts deviation coding for each level in
#' comparison with the grand mean.
#' \code{Codingtype = "standard"} conducts standard dummy coding for each level
#' in comparison with the reference level (first level).
#'
#' @param Coef_initial A \eqn{p}-dimensional vector for the initial coefficient 
#' value of the IHT algorithm.  The default is to use Lasso with the sparsity 
#' closest to \eqn{n-1}. 
#' 
#' @param penalize_mod A logical flag to indicate whether adjustment is used in
#' ranking groups of features. This augment is applicable only when
#' \code{categorical= TRUE} with \code{group=T}; the default is true:
#' a factor of \eqn{\sqrt{J}} is divided from the \eqn{L_2} effect of a group with J members.
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
#' @param max_iter Maximum number of iteration steps. Default is 500. 
#'
#' @param tol A tolerance level to stop the iteration, when the squared sum of
#' differences between two successive coefficient updates is below it.
#' Default is \eqn{10^{-2}}.
#' 
#'
#' @param selection A logical flag to indicate whether an elaborate selection
#' is to be conducted by \code{smle_select} after screening.
#'  Default is FALSE. If TRUE, the function will return a '\code{selection}' object, 
#'  see \code{smle_select()} documentation for more details.
#'
#' @param ... Additional arguments to be passed to \code{smle_select()} if \code{selection=TRUE}. 
#' See \code{smle_select()} documentation for more details. 
#' 
#' @references
#' UCLA Statistical Consulting Group. \emph{coding systems for categorical
#' variables in regression analysis}. \url{https://stats.idre.ucla.edu/spss
#' /faq/coding-systems-for-categorical-variables-in-regression-analysis-2/}.
#' Retrieved May 28, 2020.
#'
#' Xu, C. and Chen, J. (2014). The Sparse MLE for Ultrahigh-Dimensional Feature
#' Screening, \emph{Journal of the American Statistical Association}, \bold{109}(507), 1257–1269.
#'
#'
#'
#'
#'
#' @return
#' Returns a \code{'smle'} object with
#' \item{I}{A list of iteration information.
#'
#' \code{Y}: Same as input Y.
#'
#' \code{CM}: Design matrix of class \code{'matrix'} for numeric features (or \code{'data.frame'} with categorical features).
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
#' The output includes both features retained by SMLE and the features specified in \code{'keyset'}.}
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
#' \item{Usearch}{Number of times in searching a proper \eqn{u^{-1}} at each step over the IHT iterations.}
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
#' fit<-SMLE(Y=Data$Y, X=Data$X, k=9, family = "gaussian")
#' fit
#' summary(fit)
#' ## The important features we missed:
#' setdiff(Data$subset_true,fit$ID_Retained)
#' ## Check if the important featrues are retained.
#' Data$index %in% fit$ID_Retained
#' plot(fit)
#'
#'
SMLE <- function (object, ...)
  UseMethod("SMLE")

#' @rdname SMLE
#' @export
SMLE.default<-function(object=NULL, X=NULL,Y=NULL,data=NULL, k=NULL, 
                       family=c("gaussian","binomial","poisson"),
                       categorical = NULL , keyset = NULL, intercept = TRUE ,
                       group = TRUE , codingtype = NULL , Coef_initial=NULL,
                       max_iter = 500 , tol = 0.01 ,selection = F ,
                       standardize = TRUE , fast = FALSE , U_rate=0.5 ,
                       penalize_mod = TRUE,...){

  #-------Input preprocess-------

  family<-match.arg(family)

  cl<-match.call()
  cl[[1]] <- as.name("SMLE")
  
  ####Input object is a design matrix
  if(is.null(Y)&is.null(object)){stop("Response required")}
  if(is.null(X)&is.null(data)){stop("Featrue Matrix required")}
  if(is.null(Y)&!is.null(object)){Y<-object}
  if(is.null(X)&!is.null(data)){X<-data}
  ####
  if(is.null(k)){
    
    k<-floor(1/2*log(dim(X)[1])*dim(X)[1]^(1/3))
  }
  
  if(is.null(categorical)){
    
    Ci<-(1:dim(X)[2])[sapply(X,is.factor)]
    
    if( any(Ci) ){
    
        categorical = TRUE
    
        }else{
        
          categorical= FALSE
        
          X<- as.matrix(X, ncol=dim(X)[2])
        
          }
    
    }


  #-------Run Algoriathm------

  if(!is.null(Coef_initial)){

    stopifnot(length(Coef_initial)==dim(X)[2])
  }


  if(categorical == TRUE){
    
    Ci<-(1:dim(X)[2])[sapply(X,is.factor)]
  
    if(is.null( codingtype )){codingtype<-"DV"}
    
    if(group==FALSE &codingtype!="all"){stop("codingtype should be 'all' when group = FALSE")}
    
    
    S<-Standardize(X[,-Ci])
    
    Xs<-S$Xs
    
    X_mean<-S$X_mean
    
    X_sd<- S$X_sd
    
    
    for(i in 1:length(Ci)){
      
      
      if( Ci[i] == 1 ){
        
        name<-dimnames(Xs)[[2]]
        
        Xs <- data.frame( X[,1] , Xs )
        
        colnames(Xs)<-c(names(X)[1],name)

      }else if( Ci[i]>dim(Xs)[2]){
        
        name<-dimnames(Xs)[[2]]
        
        Xs<-data.frame(Xs,X[,Ci[i]])
        
        colnames(Xs)<-c(name,names(X)[Ci[i]])

      }else{ 
        
        name<- dimnames(Xs)[[2]]
        
        Xs<-data.frame(Xs[,1:Ci[i]-1],X[,Ci[i]],Xs[,Ci[i]:dim(Xs)[2]])
        
        colnames(Xs)<-c(name[1:Ci[i]-1],names(X)[Ci[i]],name[Ci[i]:(length(name))])

        }
      
    }
    
    
    


    fit<- ctg_fit(Y, Xs, k, family, categorical, keyset, Ci, max_iter, tol, fast,
                  
                  intercept, group, codingtype, penalize_mod,
                  
                  U_rate, X_mean = X_mean, X_sd=X_sd, Coef_initial, cl)

    return(fit)

  }else{
  
    X_mean = NULL;X_sd=NULL;Xs<-X

    if(standardize== TRUE){
      
      S<-Standardize(X)
      
      Xs<-S$Xs
      
      X_mean<-S$X_mean
      
      X_sd<- S$X_sd

    }
    
    fit<-SMLE_fit(Y=Y, X= Xs, k= k, family=family, keyset=keyset,
                
                intercept=intercept, max_iter =max_iter, tol =tol, fast =fast,
                
                U_rate=U_rate, X_mean = X_mean, X_sd=X_sd, Coef_initial=Coef_initial)
  }

  #------Adding method-----

  fit$call=cl

  class(fit)="smle"

  #-----Stage II ---------
  if(selection==F){

    return(fit)

    }else{

  return(smle_select(fit,...))
  }
}

#' @rdname SMLE
#' @param object an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' @export
SMLE.formula<- function (object, data, categorical=NULL,...) {

  ## keep call (of generic)
  cl <- match.call()
  cl[[1]] <- as.name("SMLE")
  
  ## model frame
  mf <- match.call(expand.dots = FALSE)
  m <- c("formula", "data")
  m <- match(m, names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(model.frame)
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf, parent.frame())
  
  ## model terms
  mt <- attr(mf, "terms")
  
  ## model matrix and response
  x <- model.frame(mt, mf)[,-1]
  y <- model.response(mf, "numeric")
  
  if(is.null(categorical)){
    
    Ci<-(1:dim(x)[2])[sapply(x,is.factor)]
    
    if( any(Ci) ){
      
      ans <- SMLE(Y=y, X=x, categorical = TRUE,...)
      
    }else{
      
      x<- as.matrix(x, ncol=dim(x)[2],categorical= FALSE,...)
      
      
      ## fit subsets
      ans <- SMLE(Y=y, X=x, ...)
    }
  }

  
  ## done
  ans
}

SMLE_fit<-function(Y , X, k, family="gaussian", keyset=NULL,
                   intercept=TRUE, max_iter=500, tol=0.01, fast=FALSE,
                   U_rate=0.5,X_mean=NULL,X_sd=NULL,Coef_initial=NULL){
  
  LH<-rep(0,max_iter)
  
  number_of_Ucheck<-rep(0,max_iter)
  
  n<-dim(X)[1];p<-dim(X)[2]
  
  if( is.null( Coef_initial ) ){
    
    fit_pre<-glmnet(x=X,y=Y,family=family)
    
    Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
    
  }else{ 
    
    Beta0 <- Coef_initial }
  
  FD<-NULL
  
  if(is.null(keyset)){
    
    number_of_ID_retained <-k
    
  }
  else{
    
    number_of_ID_retained<-k-(length(keyset)) 
    
  }
  
  if(intercept==T){
    
    Beta0<-c(mean(Y-tcrossprod(as.matrix(X),t(Beta0))),Beta0)
    
    names(Beta0)[1]="(intercept)"
    
    X_iter<-cbind(matrix(1,nrow  = n, ncol = 1),X)
    
    ID_None0<-(1:(p+1))[Beta0!=0]
    
    keyset<-c(1,keyset+1)
    
  }else{
    
    X_iter<-X
    
    ID_None0<-(1:p)[Beta0!=0]
    
  }
  
  pp<-p+intercept
  
  Screening_index<-sub_off(1:pp,keyset)
  
  beta_path<-as.matrix(Beta0,nrow=pp,ncol=1)
  
  I<-list(CM=X,Y=Y,IM=X_iter,Beta0=Beta0,family=family)
  
  pp<-p+intercept
  
  coef_None0 <- as.matrix(Beta0[ID_None0],ncol=1)
  
  Xs_0 <- X_iter[, ID_None0]
  
  if(!is.null(Coef_initial) & sum(Coef_initial!=0)==0){ 
    
    R_0<-matrix(0, ncol=1, nrow=n) 
    
  }else{
    R_0  <- Xs_0 %*% coef_None0
  }
  
  
  R_0<-switch(family,
              
              "gaussian"=R_0,
              
              "poisson"=exp(R_0),
              
              'binomial'=exp(R_0)/(1+exp(R_0)))
  
  V_0<-crossprod(X_iter,Y - R_0)
  
  
  if(!is.null(Coef_initial) & sum(Coef_initial!=0)==0){
    
    if(intercept==T){
      
      #When starts from zero with uu is 1/||X||\infit = 1/sqrt(n) 
      uu<-1/(sqrt(n))
      
      
    }else{
      
      uu<- 1/max(colSums(Xs_0^2))
      
    }
    
    
  }else{    
    
    uu<-1/max(colSums(Xs_0^2))
    
  }
  
  
  
  U_i <- uu
  
  #iteration start form 1
  i<-1
  
  Beta_s<-Beta0
  
  
  #######################################################################
  # iteration start---------------------------
  
  
  repeat{
    
    count<-0
    
    repeat{
      
      Beta_t<- Beta_s + uu * V_0
      
      Beta_t[Screening_index]<-Hard(t=Beta_t[Screening_index],k= number_of_ID_retained)
      
      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)
      
      if ( ucheck >= 1 ){
        
        break
        
      }else{
        
        uu <- U_rate * uu
        
        count<-count+1
        
      }
      
    }
    
    sindex<-(1:pp)[Beta_s!= 0]
    
    tindex<-(1:pp)[Beta_t!= 0]
    
    fs<-sum(!(tindex %in% sindex))
    
    FD<-c(FD,fs)
    
    beta_path<-cbind(beta_path,as.matrix(Beta_t,ncol=1,nrow=pp))
    
    LH[i]<-lh(Y, X_iter, Beta_t,family=family)
    
    valid_LH_diff<- 0.01 *(LH[2] - LH[1])
    
    number_of_Ucheck[i]<-count
    
    
    ######## convergence check #######################################
    if(i>1){
      
      MSE<- sqrt(sum(( Beta_s-Beta_t )^2))
      
      if(fast == TRUE){
        
        
        if( MSE/number_of_ID_retained < tol){break}
        else if( (LH[i]-LH[i-1])< valid_LH_diff){break}
        else if(i>10){if(sum(tail(FD,10))==0){break}}
        
      }else{
        
        if(MSE < tol || i >= max_iter){
          break}
      }
      
      
    }
    
    ##################################################################
    
    
    Beta_s<-Beta_t
    
    ID_None0<- (1:pp)[Beta_s!= 0]
    
    coef_None0<-as.matrix(Beta_s[ID_None0],ncol=1)
    
    Xs_0<-X_iter[, ID_None0]
    
    R_0<-crossprod(t(Xs_0),coef_None0)
    
    R_0<-switch(family,
                
                "gaussian"=R_0,
                
                "poisson"=exp(R_0),
                
                'binomial'=exp(R_0)/(1+exp(R_0))
    )
    
    V_0<-crossprod(X_iter,Y - R_0)
    
    uu<-U_i
    
    i<-i+1
    
  }
  
  
  Coef_dist<-apply((beta_path[,-1]-beta_path[,-ncol(beta_path)]),2,function(x){sqrt(sum(x^2))})
  
  # Rescale all output
  
  if( intercept == T ){
    
    #rescale beta_path
    
    beta_path<-beta_path[-1,]             
    
    Intercept_value<- as.vector(coef_None0)[1]
    
    ID_None0<-ID_None0[-1]-1
    
    coef_None0 <- as.vector(coef_None0)[-1]
    
    if(!is.null(X_mean)){
      
      coef_None0<- coef_None0/X_sd[ID_None0]
      
      Intercept_value <- Intercept_value-X_mean[ID_None0]%*%coef_None0
      
    }
  }else{
    
    if(!is.null(X_mean)){
      
      coef_None0<- coef_None0/X_sd[ID_None0]
      
      beta_path<- beta_path/X_sd
    }
    
    Intercept_value = NULL
    
  }
  
  
  
  
  fit<-list(I=I,X=I$CM, Y=I$Y,
            
            keyset = keyset, family = family, k = k,
            
            Intercept=Intercept_value,
            
            steps = i,
            
            LH=LH[1:i],
            
            Usearch=number_of_Ucheck[1:i],
            
            Path_Retained  = beta_path,
            
            Num_Retained   = length(ID_None0),
            
            ID_Retained    = ID_None0,
            
            Coef_Retained  = coef_None0,
            
            FD=FD, ctg = FALSE, Coef_dist=Coef_dist,
            
            fast=fast
  )
}

