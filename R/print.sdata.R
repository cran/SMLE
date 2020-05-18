#' Print function for simulated data
#' @description This functions prints a summary of a data set sumulated with Gen_data.
#' In particular, it prints which features were selected to be causal, their coefficients and about the correlation structure.
#'
#' @param x \code{"sdata"} object from Gen_Data function.
#' @param ... Other parameter for print function.
#' @return
#' No return value, called for side effects.
#' @export
#' @method print sdata
#' @examples
#' Data<-Gen_Data(family ="binomial",correlation = "ID")
#' cor(Data$X[,1:10])
#' print(Data)
print.sdata<-function(x,...){


  Description<-data.frame("Dim_of_Y" = paste(c(length(x$Y),1),collapse = ' x '),
                "Dim_of_X" = paste(dim(x$X),collapse = ' x '),
                "Model_Type"= .simpleCap(x$family),
                "Relevant_features"=paste("V",x$index[order(x$index, decreasing = FALSE)],sep='',collapse=' '),
                "Coefficient" = paste(format(x$Beta, digits = 3),collapse = ' '),
                "Correlation" = .simpleCap(x$correlation)
                )

  if(x$Cate==TRUE){Description<-cbind(Description,
                                      "Categorical_features"=paste("C",x$CI,sep='',collapse=' '),
                                      "Level_of_categories" =paste(x$levels,collapse=' ')
                                      )}

  print(t(Description),...)
}

