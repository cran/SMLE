#' Prediction based on SMLE screening and selection
#'
#' @description
#' Similar to the usual predict methods,this function returns predicted mean
#' values of the response based on the features retained in \code{'smle'} object
#' or selected by \code{'selection'} object.
#'
#' @importFrom stats binomial gaussian glm.fit poisson predict.glm
#'
#'
#' @param object A fitted object of class \code{'smle'} , as the output from
#' SMLE; or \code{'selection'} as the output from smle_select.
#'
#'
#' @param newdata Matrix of new values for x at which predictions are to be made,
#' without the intercept term. If omitted, the fitted linear features are used.
#'
#'
#' @param type  Type of prediction required. "response" gives fitted values for
#' "gaussian"; fitted probabilities for 'binomial', fitted mean for 'poisson'.
#' "link" returns prediction on the scale of the linear predictors. (Same to
#' "response" in "gaussian" models)
#'
#' @param ... 	Further arguments pass to predict.glm().
#'
#'
#'
#' @return
#' Returns a vector of the predicted mean values of the response based on
#' 'newdata'and the features retained in 'object'. The predicted values depend on the
#' model specified in 'type'.
#'
#' @examples
#'
#' set.seed(123.456)
#'
#' Data_sim<-Gen_Data(n= 200, p =1000, correlation="AR",family = "gaussian")
#'
#' fit<-SMLE(Data_sim$Y,Data_sim$X, family = "gaussian")
#'
#' predict(fit , type ="link")
#'
#' @export
predict.smle<-function(object, newdata = NULL,type = c("link", "response"),...){


  if( is.null(newdata) ){

    X_s<-object$I$CM[,object$ID_Retained]

    }else{X_s <-newdata[,object$ID_Retained]}

  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())

  if( any(sapply(X_s,is.factor ))  ){

    X_dummy<- as.matrix(suppressWarnings(dummy.data.frame(X_s ,sep="_")))

    gfit<-glm.fit(X_dummy,object$I$Y,family = family)

  }else{

    gfit<-glm.fit(X_s,object$I$Y,family = family)

  }

  return(predict.glm(gfit,type = type,...))


}
#' @rdname  predict.smle
#' @examples
#'
#' E<-smle_select(fit, tune="ebic")
#'
#' predict(E , type ="link")
#'
#' @export
predict.selection<-function(object,newdata = NULL,type = c("link", "response"),...){


  if( is.null(newdata)){

    X_s<-object$X[,object$ID_Selected]

  }else{X_s <-newdata[,object$ID_Selected]}

  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())

  if( any(sapply(X_s,is.factor ))  ){

    X_dummy<- as.matrix(suppressWarnings(dummy.data.frame(X_s ,sep="_")))

    gfit<-glm.fit(X_dummy,object$Y,family = family)

  }else{

    gfit<-glm.fit(X_s,object$Y,family = family)

  }

  return(predict.glm(gfit,type = type,...))


}
