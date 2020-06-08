#' Print a selection object from smle_select
#' @description This function prints a summary of a \code{'selection'} object.
#' In particular, it gives the selected features along with their re-fitted model coefficients.
#' For reference, it also shows the values of the selection criterion used in selection for all candidate models.
#'
#' @importFrom utils capture.output
#'
#' @param x Fitted \code{'selection'} object.
#'
#' @param ... Other parameter to print.
#'
#' @return
#' No return value, called for side effects.
#'
#' @export
#'
#'
#' @method print selection
#'
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' E<-smle_select(fit)
#' print(E)
#'
print.selection<-function( x,... ){

  if(x$vote != T){

    Description<-data.frame("Selection_criterion" = as.character(x$criterion),

                            "Selected_features" = paste("V",x$ID_Selected,sep='',collapse=' '))

    if(x$criterion=='ebic') {

      Description<-cbind(Description,"Gamma_for_ebic"=as.character(x$gamma_ebic))

    }

  }else{

    Description<-data.frame("Gamma_candidates"=as.character(x$gamma_seq),

                                "features_selected_by_voting"=paste("V",x$ID_Voted,sep='',collapse=' '))

  }


  message(paste0(capture.output(Description), collapse = "\n"))
}
