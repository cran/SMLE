#' Print a SMLE object from SMLE
#'
#' @description This functions prints a summary of a SMLE object.
#' In particular, it shows the features retained after SMLE-screening
#' and the related convergence information.
#'
#' @param x Fitted '\code{smle}' object.
#'
#' @param ... This argument is not used and listed for method consistency.
#'
#' @return
#' No return value, called for side effects.
#' @export

#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' print(fit)
#'
print.smle=function(x , ...){

  cat("\nCall: ", deparse(x$call), "\n\n")

  Description<-data.frame("Dim_of_Y :" = paste(c(dim(x$I$CM)[1],1),collapse = ' x '),

                          "Dim_of_X :" = paste(dim(x$I$CM),collapse = ' x '),

                          "Model_Type :"= .simpleCap(x$family),

                          "Retained_model_size :"= as.character(x$Num_Retained),

                          "Retained_features :"=paste("V",x$ID_Retained,sep='',collapse=' '),

                          "Coefficients :" = paste(format(x$Coef_Retained, digits = 3),collapse = ' '),

                          "Number_of_steps :" = as.character(x$steps))

  if( !is.null(x$Intercept) ){

    Description <- cbind ( Description, "Intercept :" = as.character(x$Intercept) )

    }

  if("ctg" %in% class(x)){

    Description <- cbind(Description,
                         "Categorical_features :" =paste("C", x$I$CI, sep='',collapse = ' '      ),
                          "Level_of_categories  :" =paste(sapply(x$I$CM[,x$I$CI],nlevels) , collapse = ', '))
    }

  message(paste0(capture.output(t(Description)), collapse = "\n"))

}

