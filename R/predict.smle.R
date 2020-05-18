predict.smle<-function(object,newx,type=c("response")){
  if(dim(newx)[2]!=dim(object$I$CM)[2]){stop("new_data should has the same number of features with old_data")}

  }
