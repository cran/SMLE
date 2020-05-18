#' Plots to visualize the selection steps
#'
#' @description
#' This function constructs a sparsity vs. selection criterion curve for a selection object.
#'  When EBIC is used with voting, it also constructs a histogram showing the voting result.
#' @param x Fitted \code{"selection"} object from smle_select.
#' @param ... Other graphical parameters to plot.
#' @method plot selection
#' @return
#' No return value, called for side effects.
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' E<-smle_select(fit)
#' #Then E is a object of "selection"
#' plot(E)
#'
#' @export
plot.selection<-function(x,...){
  dev.new()
  plot(x$Criterion_value,xlab="Model sparisty", ylab= paste(x$criterion,"value"))
  if(x$vote ==TRUE){
    dev.new()
    percent <- function(x, digits = 2, format = "f", ...) {
      paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
    }
    y<-data.frame("Proportion"= sort(summary(x$Voting_Retained_Feature_IDs),decreasing = T)/20)
    barplot(y$Proportion,names.arg = as.numeric(names(summary(x$Voting_Retained_Feature_IDs))),
            xlab = "Retrained Features IDs",ylab="Featrues Voting Proportion",main="Voting results",
            )
  }
}
