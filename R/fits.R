GPBnet<-function(Y,X,k,family,keyset,intercept,max_iter,tol,fast,U_rate){
  if(fast==TRUE){return(fit<-GPBnet_fast(Y,X,k,family,keyset,intercept,max_iter,tol,U_rate))}
  n<-dim(X)[1];p<-dim(X)[2]
  fit_pre<-glmnet(x=X,y=Y,family=family)
  Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
  FD<-NULL
  if(intercept==T){
    Beta0<-c(mean(Y-tcrossprod(as.matrix(X),t(Beta0))),Beta0)
    names(Beta0)[1]="(intercept)"
    X_iter<-cbind(matrix(1,nrow  = n, ncol = 1),X)
    I<-list(CM=X,Y=Y,IM=X_iter,Beta0=Beta0,family=family)
    #--------------------------------------------------------------#
    #--------------------------------------------------------------#
    ind_0<-(1:(p+1))[Beta0!=0]
    ind_0<-ind_0[-1]
    keyset<-c(1,keyset+1)
  }else{
    X_iter<-X
    ind_0<-(1:p)[Beta0!=0]
  }

  pp<-p+intercept
  bc_0 <- as.matrix(Beta0[ind_0],ncol=1)
  Xs_0 <- X_iter[, ind_0]
  R_0  <- tcrossprod(Xs_0, t(bc_0))
  R_0<-switch(family,
              "gaussian"=R_0,
              "poisson"=exp(R_0),
              'binomial'=exp(R_0)/(1+exp(R_0)))
  V_0<-crossprod(X_iter,Y - R_0)
  uu<-1/svd(Xs_0)$d[1]
  ###########################################################
  # iteration start---------------------------
  U_i <- uu
  i<-1
  Beta_s<-Beta0# starting iteration.
  beta_path<-matrix(0,nrow=pp,ncol=1)
  LH<-rep(0,max_iter)
  number_of_Ucheck<-rep(0,max_iter)
  Screening_index<-sub_off(1:pp,keyset)
  repeat{
    count<-0
    repeat
    {
      Beta_t<-Beta_s + uu * V_0
      Beta_t[Screening_index]<-Hard(t=Beta_t[Screening_index],k=k)
      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)
      if(ucheck >= 1)
      {break}else{
        uu <- U_rate * uu
        count<-count+1
      }
    }
    sindex<-(1:pp)[Beta_s!= 0]
    tindex<-(1:pp)[Beta_t!= 0]
    fs<-sum(!(tindex %in% sindex))
    FD<-cbind(FD,fs)
    beta_path<-cbind(beta_path,as.matrix(Beta_t,ncol=1,nrow=pp))
    LH[i]<-lh(Y, X_iter, Beta_t,family=family)
    number_of_Ucheck[i]<-count
    ######## convergence check ##############

    MSE<- sqrt(sum((Beta_s-Beta_t)^2))
    if((i>max_iter)||(MSE < tol) ) {break}
    #########################################
    Beta_s<-Beta_t
    ind_0<- (1:pp)[Beta_s!= 0]
    bc_0<-as.matrix(Beta_s[ind_0],ncol=1)
    Xs_0<-X_iter[, ind_0]
    R_0<-crossprod(t(Xs_0),bc_0)
    R_0<-switch(family,
                "gaussian"=R_0,
                "poisson"=exp(R_0),
                'binomial'=exp(R_0)/(1+exp(R_0)))
    V_0<-crossprod(X_iter,Y - R_0)
    uu<-U_i
    i<-i+1

  }

  index<-sub_off(ind_0,keyset)
  if(intercept==T){
    Intercept_value<-bc_0[1]
    index<-index-1
    ind_0<-ind_0[-1]-1
  }

  bc_1<-as.matrix(Beta_s[index+1],ncol=1)

  fit<-list(I=I,
            User_Selected_index=keyset,family=family,k=k,Intercept=Intercept_value,
            steps = i,LH=LH[1:i],Uchecks=number_of_Ucheck[1:i],Retained_Features_path=beta_path,
            Number_of_Retained_Features=length(ind_0),Retained_Feature_IDs=ind_0,
            Coefficients_of_Retained_Features=bc_0,FD=FD,fast=fast)
}


GPBnet_fast<-function(Y,X,k,family,keyset,intercept,maxit,tol,U_rate){
  n<-dim(X)[1];p<-dim(X)[2]
  fit_pre<-glmnet(x=X,y=Y,family=family)
  Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])

  if(intercept==T){
    Beta0<-c(mean(Y-tcrossprod(as.matrix(X),t(Beta0))),Beta0)
    names(Beta0)[1]="(intercept)"
    X_iter<-cbind(matrix(1,nrow  = n, ncol = 1),X)
    I<-list(CM=X,Y=Y,IM=X_iter,Beta0=Beta0,family=family)
    #--------------------------------------------------------------#
    #--------------------------------------------------------------#
    ind_0<-(1:(p+1))[Beta0!=0]
    ind_0<-ind_0[-1]
    keyset<-c(1,keyset+1)
  }else{
    X_iter<-X
    ind_0<-(1:p)[Beta0!=0]
  }

  pp<-p+intercept
  bc_0 <- as.matrix(Beta0[ind_0],ncol=1)
  Xs_0 <- X_iter[, ind_0]
  R_0  <- tcrossprod(Xs_0, t(bc_0))
  R_0<-switch(family,
              "gaussian"=R_0,
              "poisson"=exp(R_0),
              'binomial'=exp(R_0)/(1+exp(R_0)))
  V_0<-crossprod(X_iter,Y - R_0)
  uu<-1/max(rowSums(crossprod(Xs_0)))
  ###########################################################
  # iteration start---------------------------
  U_i <- uu
  i<-1
  Beta_s<-Beta0# starting iteration.
  Screening_index<-sub_off(1:pp,keyset)
  repeat{
    repeat
    {
      Beta_t<-Beta_s + uu * V_0
      Beta_t[Screening_index]<-Hard(t=Beta_t[Screening_index],k=k)
      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)
      if(ucheck >= 1)
      {break}else{
        uu <- U_rate * uu
      }
    }
    ######## convergence check ##############
    MSE<- sqrt(sum((Beta_s-Beta_t)^2))/k
    if((i>maxit)||(MSE < tol) ) {break}
    #########################################
    Beta_s<-Beta_t
    ind_0<- (1:pp)[Beta_s!= 0]
    bc_0<-as.matrix(Beta_s[ind_0],ncol=1)
    Xs_0<-X_iter[, ind_0]
    R_0<-crossprod(t(Xs_0),bc_0)
    R_0<-switch(family,
                "gaussian"=R_0,
                "poisson"=exp(R_0),
                'binomial'=exp(R_0)/(1+exp(R_0)))
    V_0<-crossprod(X_iter,Y - R_0)
    uu<-U_i
    i<-i+1

  }

  index<-sub_off(ind_0,keyset)
  if(intercept==T){
    Intercept_value<-bc_0[1]
    index<-index-1
    ind_0<-ind_0[-1]-1
  }

  bc_1<-as.matrix(Beta_s[index+1],ncol=1)

  fit<-list(I=I,keyset=keyset,family=family,k=k,Intercept=Intercept_value,steps=i,
            Number_of_Retained_Features=length(ind_0),Retained_Feature_IDs=ind_0,
            Coefficients_of_Retained_Features=bc_0,fast=TRUE)
}
