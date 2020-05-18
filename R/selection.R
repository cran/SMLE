#‘ Elaborative feature selection with SMLE
#'
#' Elaborative feature selection with SMLE
#'
#' @description Given a response and a set of \code{K} features, this function first runs \code{SMLE (fast=TRUE)} to generate a series of sub-models with sparsity \code{k} varying from \code{k_min} to \code{k_max}.
#' It then selects the best model from the series based on a selection criterion.
#' When criterion EBIC is used, users can choose to repeat the selection with different values of the tuning parameter, gamma, and conduct importance voting for each feature.
#'
#' @details
#' There are three types of input allowed:
#' Object with class "smle", the output from main function SMLE;
#' Object with class "sdata", the ouput from Gen_Data;
#' Input data pair directly by Y, X. It is not recommender to use object of type sdata or the data matrices X,Y for high demensional data.
#'
#' @references
#'
#' Chen. J. and Chen. Z. (2012). "Extended BIC for small-n-large-P sparse GLM." \emph{Statistica Sinica}: 555-574.
#'
#' Chen. J. and Chen. Z. (2008). "Extended Bayesian information criteria for model selection with large model spaces."  \emph{Biometrika} 95.3: 759-771.
#'
#' Chen, Z. and Chen. J. (2009). "Tournament screening cum EBIC for feature selection with high-dimensional feature spaces." \emph{Science in China Series A: Mathematics} 52.6 : 1327-1341.
#'
#' @importFrom parallel detectCores mclapply
#'
#' @param x Object of class "smle" or "sdata", or directly input data pair (Y,X).
#' @param ... Other parameters.
#'
#' @return
#' Returns a \code{"selection"} object with
#' \item{Retained_Feature_IDs}{A list of varible IDs selected.}
#' \item{Coeffients_of_Retained_Features}{A list of coefficients for selected features fit by glmnet }
#' \item{Criterion_value}{A list of value according to selected criteria and model sparisity.}
#' \item{Voting_Retained_Feature_IDs}{A list of Voting selection results; item returned only when vote==T}
#'
#' @examples
#'
#'# This a simple example for Gaussian assumption.
#'Data<-Gen_Data(correlation="MA",family = "gaussian")
#'fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#'E<-smle_select(fit)
#'plot(E)
#' @export
#'
smle_select<-function(x, ...){
  UseMethod("smle_select")
}



#' @rdname  smle_select
#' @method smle_select smle
#' @S3method smle_select smle
#' @export
smle_select.smle<-function(x,...){
  Data<-structure(list(Y=x$I$Y,X=x$I$CM,family=x$I$family),class = "sdata")
  S<-smle_select(Data,sub_model =x$Retained_Feature_IDs,...)
  return(S)
}


#' Elaborative feature selection with SMLE
#' @rdname smle_select
#' @method smle_select sdata
#' @S3method smle_select sdata
#' @import doParallel foreach
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#'
#' @param k_min The lower bound of target model sparsity. Default is 1.
#' @param k_max The upper bound of target model sparsity. Default is as same as the number of columns in input.
#' @param sub_model A subset of columns indicating that which columns are able to be selected.(Only for object of "sdata" and (Y,X) pair)
#' @param gamma_ebic Parameter for Extended Bayesian Information Criteria. Must be v between (0, 1). Default is 0.5.
#' @param vote The logical flog for whether to perform the voting procedure. Only available when \code{tune ='ebic'}.
#' @param tune Selection criterion, must bu one of  \code{'aic','bic', or 'ebic'}. Default is 'ebic'.
#' @param gamma_seq The sequence of values for gamma_ebic when \code{vote =TRUE}.
#' @param vote_threshold A relative voting threshold in percentage. A feature is considered to be important when it receives votes passing the threshold.
#' @param para Logical flag to use parallel computing to do voting selection. Default is FALSE. see Details.
#' @param num_cores The number of cores to use. The default will be all cores detected.
#' @export
smle_select.sdata<-function(x, k_min=1, k_max=10, sub_model=NULL,
                            gamma_ebic=0.5, vote= FALSE, tune=c('ebic','aic','bic'),
                            gamma_seq=c(seq(0,1,0.2)), vote_threshold=NULL, para = FALSE, num_cores=NULL,...){
  #input check
  tune=match.arg(tune)
  X<-x$X
  Y<-x$Y
  family<-x$family
  n<-dim(X)[1]
  pp<- dim(X)[2]
  #check input data dimension, warning for High dimensional input.
  if(is.null(sub_model)){
    if(pp>1000 & pp>n ){
      input <- readline("For high dimensional data input, we recommand to do screening(smle) before selection,\n Do you still want to continue the funciton? Stop it by typing no ")
      if (substr(input, 1, 1) == "n"){stop("Algorithm stopped by user")}
      else{if(vote==T){
        cat("Directly feature selection by voting ebic is strongly not recommended, its relatively slow \n
                           please be patient, and you may want to try parallel by setting para = TURE \n")
        input <- readline("Stop it by typing: no ")
        if (substr(input, 1, 1) == "n"){stop("Algorithm stopped by user")}
      }
      } }
    X_s<-X
  }else{X_s<-X[,sub_model]}
  if(is.null(vote_threshold)){vote_threshold=(k_min+k_max)/2}
  IP<-NULL;Voting_Retained_Feature_IDs<-NULL;vs<-NULL
  if(para==TRUE){
    if(is.null(num_cores)){
      num_cores<-parallel::detectCores()
    }
    if(.Platform$OS.type == "windows"){cat("We are using windows multi-cores computing with", num_cores,"cores. \n")}
    if(.Platform$OS.type == "unix"){cat("We are using unix multi-cores computing with", num_cores,"cores. \n")}
  }
  if(k_min<0 || k_max<k_min ||k_min%%1!=0 ||k_max%%1 !=0||k_max>pp){stop("Retained model size setting error(k_min and k_max).")}

  #Feature selection by criterion.
  criter_value<-ebicc(Y,X_s,family,tune,k_min,k_max,n,pp,gamma_ebic,para,num_cores)
  v_s <- which.min(criter_value)+k_min-1
  f_s<- SMLE(Y=Y, X=X_s, k=v_s, family=family)
  #Feature selection by ebic voting.
  ebic_selection_by_gamma<-function(gamma){
    ebic_value<-ebicc(Y,X_s,family,tune,k_min,k_max,n,pp,gamma,para,num_cores)
    v_s <- which.min(ebic_value)+k_min-1
    return(SMLE(Y=Y, X=X_s, k=v_s, family=family)$Retained_Feature_IDs)
  }
  if(vote==T){
    stopifnot(tune == 'ebic')
    vs<-c()
    if(para==TRUE){
        vs<-unlist(mclapply(gamma_seq,ebic_selection_by_gamma))
    }else{vs<-unlist(lapply(gamma_seq,ebic_selection_by_gamma))}
    IP<-as.factor(vs)
    Voting_Retained_Feature_IDs<-as.numeric(names(summary(IP)[order(summary(IP),decreasing= T)[1:min(length(summary(IP)),vote_threshold)]]))
  }
  if(is.null(sub_model)){
    Retained_Feature_IDs=f_s$Retained_Feature_IDs
  }else{
    Retained_Feature_IDs<-sub_model[f_s$Retained_Feature_IDs]
    Voting_Retained_Feature_IDs<-sub_model[Voting_Retained_Feature_IDs]
  }


  S<-list(Retained_Feature_IDs=Retained_Feature_IDs,
          Coeffients_of_Retained_Features=f_s$Coefficients_of_Retained_Features,
          vote=vote,criterion=tune,Criterion_value=criter_value,Voting_Retained_Feature_IDs=Voting_Retained_Feature_IDs,
          gamma_ebic=gamma_ebic,gamma_seq=gamma_seq)
  class(S)<-"selection"
  return(S)
  }

#' @rdname smle_select
#' @method smle_select default
#' @S3method smle_select default
#' @param X Input features matrix.
#' @param family Response type (see SMLE); default is gaussian. When input object is smle or sdata, the same model will be used in the selection step.
#' @export
smle_select.default<-function(x, X=NULL, family='gaussian',...){
  Data<-structure(list(Y=x,X=X,family=family),class = "sdata")
  S<-smle_select(Data,...)
  return(S)
}







