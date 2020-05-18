ctg_fit<-function(Y,X,k=NULL,family=c("gaussian","binomial","poisson"),
                  categorical= FALSE ,User_index = FALSE,max_iter= 500,tol=10^(-3),
                  intercept = TRUE,group= TRUE,codingtype="DV",penalize_mod=TRUE){
  family<-match.arg(family)
  call<-match.call()
  n<-dim(X)[1];p<-dim(X)[2]
  #--------------------------------------------------------------#
  Ci<-(1:p)[sapply(X,is.factor)]
  if(codingtype=="all"){
    nlevel<-sapply(X[,Ci],nlevels)
  }else{nlevel<-sapply(X[,Ci],nlevels)-1}
  X_dummy<-suppressWarnings(dummy.data.frame(X,sep="_",codingtype = codingtype))
  #--------------------------------------------------------------#
  fit_pre<- glmnet(x=as.matrix(X_dummy),y=Y,family=family)
  Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
  Beta0<-c(mean(Y-tcrossprod(as.matrix(X_dummy),t(Beta0))),Beta0)
  names(Beta0)[1]="(intercept)"
  #-----------------d---------------------------------------------
  Dummy_index<-c()
  Dummy_sum<-0
  for(i in 1:length(Ci)){
    Dummy_index<-c(Dummy_index,list(Ci[i]+seq(nlevel[i])+Dummy_sum))
    Dummy_sum<-Dummy_sum+nlevel[i]-1
  }
  DFI<-Dummy_index

  DI<-unlist(lapply(DFI, function(l) l[[1]]))
  #--------------------------------------------------------------#
  X_iter<-as.matrix(cbind(matrix(1,nrow  = n, ncol = 1),X_dummy))
  colnames(X_iter)[1]<-'(intercept)'
  pp<-dim(X_iter)[2]
  I<-list(Y=Y,CM=X,CI=Ci,nlevel=nlevel,DM=as.matrix(X_dummy),IM=X_iter,
          DFI=DFI,DI=DI,Beta0=Beta0,family=family,codingtype=codingtype)

  #--------------------------------------------------------------#


  #--------------------------------------------------------------#
  ind_0<-(1:pp)[Beta0!=0]
  ind_0<-ind_0[-1]
  bc_0 <- as.matrix(Beta0[ind_0],ncol=1)
  Xs_0 <- X_iter[, ind_0]
  R_0  <- tcrossprod(Xs_0, t(bc_0))
  R_0<-switch(family,
              "gaussian"=R_0,
              "poisson"=exp(R_0),
              'binomial'=exp(R_0)/(1+exp(R_0))
  )
  V_0<-crossprod(X_iter,Y - R_0)
  uu <-1/max(rowSums(crossprod(Xs_0)))
  ###########################################################
  # iteration start---------------------------
  U_i<- uu
  i<-1
  Beta_s<-Beta0# starting iteration.
  beta_path<-matrix(0,nrow=pp,ncol=1)
  LH<-list()
  number_of_Ucheck<-list()
  Screening_index<-sub_off(1:p,User_index)
  Screening_Dindex<-sub_off(1:pp,CI2DI(I,User_index))[-1]
  repeat{
    count<-0
    repeat
    {
      Beta_t<-Beta_s + uu * V_0
      if(group==T){
        Beta_t<-GroupHard(Beta_t,I,k=k,Screening_index,penalize_mod)
      }else{
        Beta_t[Screening_Dindex]<- Hard(t=Beta_t[Screening_Dindex], k=k)
      }
      ########## u-check #################
      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)
      if(ucheck >= 1)
      {break}else{
        uu <- 0.5 * uu
      }
    }

    beta_path<-cbind(beta_path,as.matrix(Beta_t,ncol=1,nrow=pp))
    likelihood<-lh(Y, X_iter, Beta_t,family=family)
    LH<-c(LH,likelihood)
    ######## convergence check ##############

    MSE<- sqrt(sum((Beta_s-Beta_t)^2))
    if( (i>max_iter) ||  (MSE < tol)) {break}

    #########################################
    Beta_s<-Beta_t
    ind_0<- (1:pp)[Beta_s!= 0]
    bc_0<-as.matrix(Beta_s[ind_0],ncol=1)
    Xs_0<-X_iter[, ind_0]
    R_0<-crossprod(t(Xs_0),bc_0)
    R_0<-switch(family,
                "gaussian"=R_0,
                "poisson"=exp(R_0),
                'binomial'=exp(R_0)/(1+exp(R_0))
    )
    V_0<-crossprod(X_iter,Y - R_0)
    uu<-U_i
    i<-i+1
  }
  Intercept_value<-bc_0[1]
  bc_0<-bc_0[-1]
  ind_0<-DI2CI(I,ind_0[-1])
  index<-sub_off(ind_0,User_index)
  bc_1<-as.matrix(Beta_t[index],ncol=1)
  bc_0<-as.matrix(Beta_t[ind_0],ncol=1)

  fit<-list(I=I,
            Retained_Feature_IDs=index,
            Coefficients_of_Retained_Features= bc_1,
            User_Selected_index=User_index,
            family=family,
            k=k,
            Intercept=Intercept_value,
            steps = i,
            LH=LH[1:i],
            Retained_Features_path=beta_path,
            Number_of_Retained_Features=length(ind_0),
            All_Nonzeros_IDs=ind_0,
            All_Nonzeroes_coefficients=bc_0,
            group=group)

  fit$call=call
  class(fit)=c(class(fit),"ctg","smle")
  fit

}
