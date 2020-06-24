ctg_fit<-function(Y , X , k ,

                  family = c("gaussian","binomial","poisson"),

                  categorical  , keyset ,

                  max_iter, tol ,

                  intercept , group ,

                  codingtype , penalize_mod,

                  call, U_rate  ){

  family<-match.arg(family)

  n<-dim(X)[1];p<-dim(X)[2]

  #--------------------------------------------------------------#

  Ci<-(1:p)[sapply(X,is.factor)]

  if(codingtype=="all"){

    dum_col<-sapply(X[,Ci],nlevels)

  }else{

    dum_col<-sapply(X[,Ci],nlevels)-1

  }

  X_dummy<-suppressWarnings(dummy.data.frame(X,sep="_",codingtype = codingtype))

  #--------------------------------------------------------------#


  fit_pre<- glmnet(x=as.matrix(X_dummy),y=Y,family=family)

  Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])

  Beta0<-c(mean(Y-tcrossprod(as.matrix(X_dummy),t(Beta0))),Beta0)

  names(Beta0)[1]="(intercept)"

  #--------------------------------------------------------------

  Dummy_index<-c()

  Dummy_sum<-0

  for(i in 1:length(Ci)){

    Dummy_index<-c(Dummy_index,list(Ci[i]+seq(dum_col[i])+Dummy_sum))

    Dummy_sum<-Dummy_sum+dum_col[i]-1

    }

  DFI<-Dummy_index

  DI<-unlist(lapply(DFI, function(l) l[[1]]))

  #--------------------------------------------------------------#

  X_iter<-as.matrix(cbind(matrix(1,nrow  = n, ncol = 1),X_dummy))

  colnames(X_iter)[1]<-'(intercept)'

  pp<-dim(X_iter)[2]

  I<-list(Y=Y,CM=X,CI=Ci,dum_col=dum_col,IM=X_iter,
          DFI=DFI,DI=DI,family=family,codingtype=codingtype)


  #--------------------------------------------------------------#




  #--------------------------------------------------------------#

  ID_None0<-(1:pp)[Beta0!=0]  #length : = pp

  coef_None0 <- as.matrix(Beta0[ID_None0],ncol=1)

  Xs_0 <- X_iter[, ID_None0]

  R_0  <- tcrossprod(Xs_0, t(coef_None0))

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

  number_of_Ucheck<-rep(0,max_iter)

  FD<-NULL

  Screening_index<-sub_off(1:p,keyset)

  Screening_Dindex<-sub_off(1:pp,CI2DI(I,keyset))

  repeat{

      count<-0

      repeat
    {


      Beta_t<-Beta_s + uu * V_0         # length(Beta_s)  = pp

      if(group==T){

          Beta_t<-GroupHard(Beta_t,I,k=k-(length(keyset)),Screening_index,penalize_mod)

          # length(Beta_t) = pp

          }else{

            Beta_t[Screening_Dindex]<- Hard(t=Beta_t[Screening_Dindex], k=k-(length(keyset)))

            }


      ########## u-check #################

      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)

      if(ucheck >= 1)
      {
        break

        }else{

          uu <- U_rate * uu

          count<-count+1

          }
    }

    sindex<-(1:pp)[Beta_s!= 0]

    tindex<-(1:pp)[Beta_t!= 0]

    fs<-sum(!(tindex %in% sindex))

    FD<-cbind(FD,fs)

    beta_path<-cbind(beta_path,as.matrix(Beta_t,ncol=1,nrow=pp))

    likelihood<-lh(Y, X_iter, Beta_t,family=family)

    LH<-c(LH,likelihood)

    number_of_Ucheck[i]<-count

    ######## convergence check ##############

    MSE<- sqrt(sum((Beta_s-Beta_t)^2))

    if( (i>max_iter) ||  (MSE < tol)) {break}

    #########################################

    Beta_s<-Beta_t

    ID_None0<- (1:pp)[Beta_s!= 0]

    coef_None0 <- as.matrix(Beta_s[ID_None0],ncol=1)

    Xs_0<-X_iter[, ID_None0]

    R_0<-crossprod(t(Xs_0), coef_None0)

    R_0<-switch(family,
                "gaussian"=R_0,
                "poisson"=exp(R_0),
                'binomial'=exp(R_0)/(1+exp(R_0))
    )

    V_0<-crossprod(X_iter,Y - R_0)

    uu<-U_i

    i<-i+1
  }



  Intercept_value<-coef_None0[1]



  if(group == TRUE){

    ID_None0<-DI2CI(I,ID_None0)   # issue here

    coef <- Group_Beta(Beta_s,I,penalize_mod)

    coef_None0 <- coef[coef!=0]

  }else{

    ID_None0<-DI2CI(I,ID_None0)

    coef <- Group_Beta(Beta_s,I,penalize_mod)

    coef_None0 <- coef[coef!=0]

  }


  fit<-list(I=I,
            ID_Retained=ID_None0,
            Coef_Retained= coef_None0,
            keyset=keyset,
            Uchecks=number_of_Ucheck[1:i],
            family=family,
            k=k,
            Intercept=Intercept_value,
            steps = i,
            LH=LH[1:i],
            Path_Retained=beta_path,
            Num_Retained   = length(ID_None0),
            group=group,
            fast =FALSE,
            FD= FD,
            Cate =TRUE,
            CT = codingtype
            )

  fit$call=call
  class(fit)="smle"
  fit

}
