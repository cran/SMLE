#' Data simulator for high-dimensional GLMs
#'
#' This function generates synthetic datasets from GLMs with a user-specified correlation structure.
#' It permits both numerical and categorical features, whose quantity can be larger than the sample size.
#'
#' @details
#'
#' Simulated data \eqn{(y_i , x_i)} for \eqn{i = 1, . . . , n} are generated as follows:
#' First, randomly sample \code{num_truecoef} important features among the \code{p} features, the magnitude of the effects from a U(0,1) and the sign of the effect randomly(\eqn{\beta}'s).
#' Second, generate X using theselected correlation structure (independent,  auto-regressive, moving average, compound symmetry). Then, to generating categorical data, we convert numerical
#' features to categorical by binning value. We randomly select num_ctgixdx numerical columns and convert them into a four level factor.
#'
#' Moving average: candidate features \eqn{x_1,..., x_p} are joint normal,
#' marginally N(0, 1), with \eqn{cov(x_j, x_{j-1}) = \frac{2}{3}}, \eqn{cov(x_j, x_{j-2}) = \frac{1}{3}} and \eqn{cov(x_j, x_h) = 0} for \eqn{|j-h| \geq 3}.
#'
#' Compound symmetry: candidate features \eqn{x_1,..., x_p} are joint normal, marginally N(0, 1), with \eqn{cov(x_j, x_h) =0.15} if \eqn{\ j\ ,\ h\ } are both in the set of important features and \eqn{cov(x_j, x_h) = 0.3} when only
#' one of \eqn{j} or \eqn{h} are in the set of important features.
#'
#' Auto-regressive: candidate features \eqn{x_1,..., x_p} are joint normal marginally N(0, 1), with \eqn{cov(x_j, x_{j+1}) = \rho} for all \eqn{j}.
#'
#' Then, generate the response variable Y according to its response type. For Gaussian model, \eqn{Y =x^T \cdot \beta + \epsilon} where \eqn{\epsilon \in N(0,1)}.
#' For the binary model let \eqn{\pi = P(Y = 1|x)}. Sample y from Bernoulli(\eqn{\pi}) where \eqn{logit(\pi) = x^T \cdot\beta}.
#' Finally, for the Poisson model, Y is generated from Poisson distribution with the link \eqn{\pi =exp(x^T \cdot \beta )}. For more details (see reference below)
#'
#' @param n Sample size, number of rows of the data.frame (or matrix).
#' @param p Number of features.
#' @param num_ctgidx The number of features that are categorical. Set to FALSE for only numerical features. Default is FALSE.
#' @param num_truecoef The number of features that affect response. Default is 5.
#' @param family Response type.
#' \code{'gaussian'} for normally distributed data, \code{'poisson'} for non-negative counts,
#' \code{'binomial'} for binary (0-1).
#' @param correlation correlation structure between features. \code{correlation = 'ID'} for all variables independent,
#' \code{correlation = 'MA'} for moving average, \code{correlation = 'CS'} for compound symmetry, \code{correlation = 'AR'} for auto correlation. Default is "independent".
#' For more information see details.
#' @param rho Parameter for AR(1) data, when correlation = "AR". Default is 0.5.
#' @param pos_ctgidx Vector of indices denoting which columns are categorical.
#' @param pos_truecoef Vector of indices denoting which columns affect the response variable.
#'
#' @references
#' Chen Xu and Jiahua Chen. (2014),
#' The Sparse MLE for Ultrahigh-Dimensional Feature Screening
#' * Journal of the American Statistical Association*109:507,
#' pages:1257-1269
#' @return
#' Returns a \code{"sdata"} object with
#' \item{Y}{Response variable vector of length \eqn{n}}
#' \item{X}{Feature matrix or Dataframe (Matrix if \code{num_ctgidx =FALSE} and dataframe otherwise)}
#' \item{index}{Vector of columns indices of X for the features that affect the response variables (causal features).}
#' \item{Beta}{Vector of effects for the causal features.}
#'
#' @export
#'
#' @examples
#' #Simulating data with binomial response and independent strcture.
#' Data<-Gen_Data(family ="binomial",correlation = "ID")
#' cor(Data$X[,1:5])
#' print(Data)
#'
#'
Gen_Data<-function(n=200,p=5000,num_ctgidx=NULL,pos_ctgidx=NULL,num_truecoef=NULL,pos_truecoef=NULL,
                   correlation=c("ID","AR","MA","CS"),rho=0.5,family=c("gaussian","binomial","poisson")){
  correlation=match.arg(correlation)
  family=match.arg(family)

  #putting default value to parameter
  if(is.null(num_ctgidx)){num_ctgidx<-0}
  if(is.null(num_truecoef)){num_truecoef<-5}
  if(is.null(pos_truecoef)){index <- sample(1:p,size=num_truecoef,replace = F)}else(index<-pos_truecoef)
  if(is.null(pos_ctgidx)){categorical_positon <- sample(1:p,size=num_ctgidx,replace = F)}else(categorical_positon<-pos_ctgidx)

  #Input parameter check
  if(num_truecoef%%1!=0){stop("Number of ground truth important features should be integer")}
  if(num_ctgidx>p){stop("The number of categorical features should be less than number of columns p")}
  if(!is.null(pos_ctgidx) & !is.null(num_ctgidx)){stopifnot(length(pos_ctgidx)==num_ctgidx)}
  if(!is.null(pos_truecoef) & !is.null(num_truecoef)){stopifnot(length(pos_truecoef)==num_truecoef)}


  sample_name<-list()
  Gene<-c("A","T","G","C")
  #Create numerical data matrix X
  numeric_data<-Gen_DesignMatrix(correlation,n,p,rho,index=index)
  if(num_ctgidx==FALSE){
    #Correlated Simulation numerical Data
    U = sample(c(-1,1),size=num_truecoef,replace= T,prob=c(0.2,0.8))
    Beta_grdtrus = (log(n)/sqrt(n)+abs(rnorm(num_truecoef))/num_truecoef)*U
    BETA<-matrix(0, nrow=p, ncol=1)
    BETA[index, 1] <- Beta_grdtrus
    theta<-tcrossprod(as.matrix(numeric_data),t(BETA))
    if(family=="gaussian"){
      Y<- theta + rnorm(n, mean=0, sd=1)
    }else if(family=="binomial"){
      pi <- exp(theta) / (1 + exp(theta))
      Y  <- rbinom(n, size=1, prob=pi)
      }else{
        mmu <- exp(theta)
        Y  <- rpois(n, lambda=mmu)
      }
    correlation=switch(correlation,
                    'ID'='independent',
                    'MA'='moving average',
                    "AR"="auto correlation",
                    'CS'='compound symmetry')
    D<-list(Y=Y,X=numeric_data,index=index,Beta=Beta_grdtrus,family=family,Cate=FALSE,correlation=correlation)
    class(D)<-"sdata"
    return(D)
    }else{
      #Converting numrical data into categorical data.
      X<-data.frame(intercept=1)
      categorical_positon<-sample(2:p,num_ctgidx,replace=F)
      categorical<-categorical_positon[order(categorical_positon)]
      Numeric2Cate<-function(data){
        Q<-quantile(data,c(0.25,0.5,0.75,1))
        for(i in 1:length(data)){
          if(as.numeric(data[i]) <Q[1]){
            data[i]<-"A"
          }else if(as.numeric(data[i])<=Q[2]){
            data[i]<-"T"
          }else if(as.numeric(data[i])<=Q[3]){
            data[i]<-"C"
          }else{data[i]<-"G"}
        }
        data
      }
      for(i in 1:p){
        if(i %in% categorical){
          sample_name[i] <- paste('Gene_pairs',i,sep="_")
          cat_data <- Numeric2Cate(numeric_data[,i])
          X<-cbind(X,cat_data)
        }else{
          sample_name[i] <- paste('numeric',i,sep="_")
          X<-cbind(X,numeric_data[,i])
        }
      }
      X<-data.frame(X)[,-1]
      colnames(X)<-sample_name
      X_dummy<-suppressWarnings(dummy.data.frame(X,sep="_",codingtype = "all"))
      pp<-dim(X_dummy)[2]

      U = sample(c(-1,1),size=length(pos_ctgidx),replace= T,prob=c(0.2,0.8))
      Beta_grdtrus = (log(n)/sqrt(n)+abs(rnorm(length(pos_ctgidx)))/length(pos_ctgidx))*U

      BETA<-matrix(0, nrow=dim(X_dummy)[2], ncol=1)
      BETA[pos_ctgidx, 1] <- Beta_grdtrus
      theta<-tcrossprod(as.matrix(X_dummy),t(BETA))
      if(family=="gaussian"){
        Y<- theta + rnorm(n, mean=0, sd=1)
        }else if(family=="binomial"){
        pi <- exp(theta) / (1 + exp(theta))
        Y  <- rbinom(n, size=1, prob=pi)
      }else{
        mmu <- exp(theta)
        Y  <- rpois(n, lambda=mmu)
        }
      #record groudtruth beta and its position
      Ci<-(1:dim(X)[2])[sapply(X,is.factor)]
      nlevel<-sapply(X[,Ci],nlevels)
      Dummy_index<-c()
      Dummy_sum<-0
      for(i in 1:length(Ci)){
        Dummy_index<-c(Dummy_index,list(Ci[i]+seq(nlevel[i])+Dummy_sum))
        Dummy_sum<-Dummy_sum+nlevel[i]-1
     }
      DFI<-Dummy_index
      DI<-unlist(lapply(DFI, function(l) l[[1]]))

      correlation=switch(correlation,
                      'ID'='independent',
                      'MA'='moving average',
                      'AR'='auto correlation',
                      'CS'='compound symmetry')
      D<-list(Y=Y,X=X,index=index,Beta=Beta_grdtrus,CI=Ci,
              family=family,Cate=TRUE,correlation=correlation,levels=nlevel)
      class(D)<-"sdata"
      return(D)
    }
}






