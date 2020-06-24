#' Plots to visualize the SMLE screening step
#'
#' This function returns two plot windows. By default, the first contains 4 plots to assess:
#' 1) log-likelihood,  2) Euclidean distance between the current
#' and the previous coefficient estimates,   3)  the number of tries in tuning parameter "u" in IHT
#' algorithm (see "Ucheck" in SMLE),  and 4) the number of features changed in the current active set.
#' By default, the second plot shows the solution path (estimated coefficient by iteration step) for selected features.
#'
#' @param x Fitted  \code{"smle"} object from SMLE.
#' @param Display For the solution path plot, show path for the most significant coefficients(\code{top_row}) or for all coefficients(\code{all}).
#' @param num_path Number of top coefficients to be shown in solution path plot if \code{type = "top_row"}. Default in solution path plot is 5.
#' @param which_path A vector to control which features are shown in addition to the paths for the most significant coefficients if \code{type ="top_row"}.
#' @param out_plot A number from 1 to 5 indicating which plot is to be shown in the separate window; the default for solution path plot is "5".
#' See Description for plot labels 1-4.
#' @param ... Additional arguments to the plot function.
#' @return
#' No return value, called for side effects.
#' @export
#' @method plot smle
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' plot(fit)
#'

plot.smle<-function(x,Display=c("top_row","all"),num_path=NULL,which_path=NULL,out_plot=5,...){


  Display <- match.arg(Display)

  oldpar <- par(no.readonly = TRUE)

  on.exit(par(oldpar))

  #-------------------check
  N_steps<-x$step

  l2_norm<-function(x){sqrt(sum(x^2))}

  l1_norm<-function(x){sum(abs(x))}

  Feature_path<-x$Path_Retained

  if(x$k>30 & Display!="top_row"){

    stop("too much features to plot, try with Display=top_row")

  }
  #--------------------plot-------------
  if(out_plot==5){

  par(mfrow=c(2,2))

  plot(y=x$LH,x=1:N_steps,xlab ="steps",ylab = 'Likelihood')

  title("Likelihood vs steps")

  plot(apply((x$Path_Retained[,-1]-x$Path_Retained[,-(N_steps+1)]),2,l2_norm)[-1],
       x=1:(N_steps-1),xlab = "steps",ylab="||\\beta_{t-1}  - \\beta_{t} ||",
       main="Euclidean distance between coefficient estimates" )

  plot(y=x$Uchecks,x=1:N_steps,xlab = "steps",ylab="number_of_Ucheck")

  title("Number of Uchecks.")

  plot(y=x$FD,x=1:N_steps,xlab ="steps",ylab="Retained feature change")

  title("Retained feature change")

  dev.new()

  if(Display =="top_row"){

    if(is.null(num_path)){num_path=5}

    num_path <-min(length(x$ID_Retained),num_path)

    TOP_value<-rep(0,5)

    TOP_value[1:num_path]<-sort(abs(Feature_path[,N_steps+1]),decreasing = T)[1:num_path]

    TOP_index<-(1:dim(x$I$CM)[2])[Feature_path[,N_steps+1]%in% cbind(TOP_value,-TOP_value)]

    TOP_index<-unique(c(TOP_index,which_path))

    Feature_path<-Feature_path[TOP_index,]

    plot(NULL, xlim=c(0,N_steps+1), ylim=c(floor(min(Feature_path[,])),ceiling(max(Feature_path[,]))),
         ylab="Coefficient", xlab="steps",...)

    title("Retained Features path")

    lines(rep(0,N_steps),lty=1,lwd=1,col="black")

    for(i in 1:length(TOP_index)){

        beta<-rep(0,N_steps)

        for(j in 1:N_steps){

          beta[j]<-Feature_path[i,j]

          }

        lines(beta,lty=i,col=rainbow(length(TOP_index))[i])
    }

    legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
           legend=TOP_index-1,bty="n")
  }else{

      plot(NULL, xlim=c(0,N_steps), ylim=c(floor(min(Feature_path[,N_steps])),ceiling(max(Feature_path[,N_steps]))),

                ylab="Coefficient", xlab="steps",...)

      title("Retained_Features_path")

      for(i in x$ID_Retained){

        beta<-rep(0,N_steps)

        for(j in 1:N_steps){

          beta[j]<-Feature_path[i,j]

          }

        lines(beta)

        }
  }
  }else if(out_plot==1){

    par(mfrow=c(2,2))
    if(Display =="top_row"){
      if(is.null(num_path)){num_path=5}
      num_path <-min(length(x$ID_Retained),num_path)
      TOP_value<-rep(0,5)
      TOP_value[1:num_path]<-sort(abs(Feature_path[,N_steps+1]),decreasing = T)[1:num_path]
      TOP_index<-(1:dim(x$I$CM)[2])[Feature_path[,N_steps+1]%in% cbind(TOP_value,-TOP_value)]
      TOP_index<-unique(c(TOP_index,which_path))
      Feature_path<-Feature_path[TOP_index,]
      plot(NULL, xlim=c(0,N_steps+1), ylim=c(floor(min(Feature_path[,])),ceiling(max(Feature_path[,]))),
           ylab="Coefficient", xlab="steps")
      title("Retained Features path")
      lines(rep(0,N_steps),lty=1,lwd=1,col="black")
      for(i in 1:length(TOP_index)){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta,lty=i,col=rainbow(length(TOP_index))[i])
      }
      legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
             legend=TOP_index-1,bty="n")
    }else{
      plot(NULL, xlim=c(0,N_steps), ylim=c(floor(min(Feature_path[,N_steps])),ceiling(max(Feature_path[,N_steps]))),
           ylab="Coefficient", xlab="steps",...)
      title("Retained_Features_path")
      for(i in x$ID_Retained){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta)
      }
    }

    plot(apply((x$Path_Retained[,-1]-x$Path_Retained[,-(N_steps+1)]),2,l2_norm),x=1:N_steps,xlab = "steps",ylab="square_difference")
    title('||beta_{t-1}  - \\beta_{t} ||_2')
    plot(y=x$Uchecks,x=1:N_steps,xlab = "steps",ylab="number_of_Ucheck")
    title("Number of Uchecks.")
    #plot(apply(x$Path_Retained[,-1],2,l2_norm),x=1:N_steps,xlab = "steps",ylab="L_2 norm")
    #title("L_2 norm of Retained Features ")
    plot(y=x$FD,x=1:N_steps,xlab ="steps",ylab="Retained feature change")
    title("Retained feature change")

  dev.new()
  plot(y=x$LH,x=1:N_steps,xlab ="steps",ylab = 'Likelihood',...)
  title("Likelihood vs steps")
  }else if(out_plot==2){
    par(mfrow=c(2,2))
    if(Display =="top_row"){
      if(is.null(num_path)){num_path=5}
      num_path <-min(length(x$ID_Retained),num_path)
      TOP_value<-rep(0,5)
      TOP_value[1:num_path]<-sort(abs(Feature_path[,N_steps+1]),decreasing = T)[1:num_path]
      TOP_index<-(1:dim(x$I$CM)[2])[Feature_path[,N_steps+1]%in% cbind(TOP_value,-TOP_value)]
      TOP_index<-unique(c(TOP_index,which_path))
      Feature_path<-Feature_path[TOP_index,]
      plot(NULL, xlim=c(0,N_steps+1), ylim=c(floor(min(Feature_path[,])),ceiling(max(Feature_path[,]))),
           ylab="Coefficient", xlab="steps",...)
      title("Retained Features path")
      lines(rep(0,N_steps),lty=1,lwd=1,col="black")
      for(i in 1:length(TOP_index)){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta,lty=i,col=rainbow(length(TOP_index))[i])
      }
      legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
             legend=TOP_index-1,bty="n")
    }else{
      plot(NULL, xlim=c(0,N_steps), ylim=c(floor(min(Feature_path[,N_steps])),ceiling(max(Feature_path[,N_steps]))),
           ylab="Coefficient", xlab="steps")
      title("Retained_Features_path")
      for(i in x$ID_Retained){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta)
      }
    }
    plot(y=x$LH,x=1:N_steps,xlab ="steps",ylab = 'Likelihood')
    title("Likelihood vs steps")
    plot(y=x$Uchecks,x=1:N_steps,xlab = "steps",ylab="number_of_Ucheck")
    title("Number of Uchecks.")
    #plot(apply(x$Path_Retained[,-1],2,l2_norm),x=1:N_steps,xlab = "steps",ylab="L_2 norm")
    #title("L_2 norm of Retained Features ")
    plot(y=x$FD,x=1:N_steps,xlab ="steps",ylab="Retained feature change")
    title("Retained feature change")

    dev.new()
    plot(apply((x$Path_Retained[,-1]-x$Path_Retained[,-(N_steps+1)]),2,l2_norm),
         x=1:N_steps,xlab = "steps",ylab="square_difference",...)
    title('||beta_{t-1}  - \\beta_{t} ||_2')
  }else if(out_plot==3){
    par(mfrow=c(2,2))
    if(Display =="top_row"){
      if(is.null(num_path)){num_path=5}
      num_path <-min(length(x$ID_Retained),num_path)
      TOP_value<-rep(0,5)
      TOP_value[1:num_path]<-sort(abs(Feature_path[,N_steps+1]),decreasing = T)[1:num_path]
      TOP_index<-(1:dim(x$I$CM)[2])[Feature_path[,N_steps+1]%in% cbind(TOP_value,-TOP_value)]
      TOP_index<-unique(c(TOP_index,which_path))
      Feature_path<-Feature_path[TOP_index,]
      plot(NULL, xlim=c(0,N_steps+1), ylim=c(floor(min(Feature_path[,])),ceiling(max(Feature_path[,]))),
           ylab="Coefficient", xlab="steps",...)
      title("Retained Features path")
      lines(rep(0,N_steps),lty=1,lwd=1,col="black")
      for(i in 1:length(TOP_index)){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta,lty=i,col=rainbow(length(TOP_index))[i])
      }
      legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
             legend=TOP_index-1,bty="n")
    }else{
      plot(NULL, xlim=c(0,N_steps), ylim=c(floor(min(Feature_path[,N_steps])),ceiling(max(Feature_path[,N_steps]))),
           ylab="Coefficient", xlab="steps")
      title("Retained_Features_path")
      for(i in x$ID_Retained){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta)
      }
    }
    plot(y=x$LH,x=1:N_steps,xlab ="steps",ylab = 'Likelihood')
    title("Likelihood vs steps")
    plot(apply((x$Path_Retained[,-1]-x$Path_Retained[,-(N_steps+1)]),2,l2_norm),x=1:N_steps,xlab = "steps",ylab="square_difference")
    title('||beta_{t-1}  - \\beta_{t} ||_2')
    #plot(apply(x$Path_Retained[,-1],2,l2_norm),x=1:N_steps,xlab = "steps",ylab="L_2 norm")
    #title("L_2 norm of Retained Features ")
    plot(y=x$FD,x=1:N_steps,xlab ="steps",ylab="Retained feature change")
    title("Retained feature change")
    dev.new()
    plot(y=x$Uchecks,x=1:N_steps,xlab = "steps",ylab="number_of_Ucheck",...)
    title("Number of Uchecks.")
  }else if(out_plot==4){
    par(mfrow=c(2,2))
    if(Display =="top_row"){
      if(is.null(num_path)){num_path=5}
      num_path <-min(length(x$ID_Retained),num_path)
      TOP_value<-rep(0,5)
      TOP_value[1:num_path]<-sort(abs(Feature_path[,N_steps+1]),decreasing = T)[1:num_path]
      TOP_index<-(1:dim(x$I$CM)[2])[Feature_path[,N_steps+1]%in% cbind(TOP_value,-TOP_value)]
      TOP_index<-unique(c(TOP_index,which_path))
      Feature_path<-Feature_path[TOP_index,]
      plot(NULL, xlim=c(0,N_steps+1), ylim=c(floor(min(Feature_path[,])),ceiling(max(Feature_path[,]))),
           ylab="Coefficient", xlab="steps")
      title("Retained Features path")
      lines(rep(0,N_steps),lty=1,lwd=1,col="black")
      for(i in 1:length(TOP_index)){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta,lty=i,col=rainbow(length(TOP_index))[i])
      }
      legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
             legend=TOP_index-1,bty="n")
    }else{
      plot(NULL, xlim=c(0,N_steps), ylim=c(floor(min(Feature_path[,N_steps])),ceiling(max(Feature_path[,N_steps]))),
           ylab="Coefficient", xlab="steps",...)
      title("Retained_Features_path")
      for(i in x$ID_Retained){
        beta<-rep(0,N_steps)
        for(j in 1:N_steps){
          beta[j]<-Feature_path[i,j]
        }
        lines(beta)
      }
    }
    plot(y=x$LH,x=1:N_steps,xlab ="steps",ylab = 'Likelihood')
    title("Likelihood vs steps")
    plot(apply((x$Path_Retained[,-1]-x$Path_Retained[,-(N_steps+1)]),2,l2_norm),x=1:N_steps,xlab = "steps",ylab="square_difference")
    title('||beta_{t-1}  - \\beta_{t} ||_2')
    plot(y=x$Uchecks,x=1:N_steps,xlab = "steps",ylab="number_of_Ucheck")
    title("Number of Uchecks.")
    dev.new()
    #plot(apply(x$Path_Retained[,-1],2,l2_norm),x=1:N_steps,xlab = "steps",ylab="L_2 norm",...)
    #title("L_2 norm of Retained Features ")
    plot(y=x$FD,x=1:N_steps,xlab ="steps",ylab="Retained feature change")
    title("Retained feature change")
  }
}

