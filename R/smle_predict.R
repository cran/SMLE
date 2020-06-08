#' Predict method for SMLE screening and selection
#'
#' @description
#' This function returns predictions based on the feature retained from
#' screening or selected by selection.
#'
#' @importFrom stats binomial gaussian glm.fit poisson predict.glm
#'
#' @details
#' First, smle_predict refits the regression model by glm.fit between response
#' variable and features selected. Then, prediction by predict.glm function.
#'
#' @param x Objects of \code{'smle'} object, as the output from SMLE; or
#' \code{'selection'} object, as the output from smle_select.
#' @param ... Other parameters.
#' @return
#' The predictor variable
#'
#' @examples
#' set.seed(123.456)
#'
#' Data_sim<-Gen_Data(n= 200, p =1000, correlation="AR",family = "gaussian")
#'
#' fit<-SMLE(Data_sim$Y,Data_sim$X, family = "gaussian")
#'
#' E<-smle_select(fit)
#'
#' p<-smle_predict(E,type = "link")
#'
#' @export



smle_predict<-function(x, ...){
  UseMethod("smle_predict")
}




#' @rdname  smle_predict
#' @method smle_predict selection
#' @param type e type of prediction required.
#' The default is on the scale of the linear predictors; the alternative
#' "response" is on the scale of the response variable. Thus for a default
#' binomial model the default predictions are of log-odds (probabilities on
#' logit scale) and type = "response" gives the predicted probabilities.
#' The "terms" option returns a matrix giving the fitted values of each
#' term in the model formula on the linear predictor scale.
#' @param terms with type = "terms" by default all terms are returned.
#' A character vector specifies which terms are to be returned
#' @export

smle_predict.selection<-function(x,type = c("link", "response", "terms"),terms = NULL,...){

  family<-switch(x$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())

  gfit<-glm.fit(x$X[,x$ID_Selected],x$Y,family = family,...)

  P<-predict.glm(gfit,type = type)

  return(P)
}



#' @rdname  smle_predict
#' @method smle_predict smle
#' @export

smle_predict.smle<-function(x, type = c("link", "response", "terms"),terms = NULL,...){

  family<-switch(x$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())

  gfit<-glm.fit(x$I$CM[,x$ID_Retained],x$I$Y,family = family,...)

  P<-predict.glm(gfit,type = type)

  return(P)
}
#'
#' @rdname smle_predict
#'
#' @method smle_predict default
#' @export
smle_predict.default<-function(x,...){

  stop("Only smle or selection object allows")


}
