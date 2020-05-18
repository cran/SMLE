.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
Uh<-function(A, uh, b0, b1, family="gaussian")
{
  if(family != "poisson")
  {
    if(family=="gaussian")
    { rho<-1 }

    if(family=="binomial")
    { rho<- 0.25 }

    bb<-b1-b0
    in_b<- (1:length(bb))[bb != 0]

    bt<-0

    if(length(in_b) != 0)
    {
      bt<- bb[in_b]
      bt<- as.matrix(bt, ncol=1)
    }

    tm1 <-  sum(bt^2) / uh
    tm2 <-  sum(tcrossprod(A[,in_b],t(bt))^2 ) * rho

    Bt <- exp(tm1 - tm2)
  }

  if(family=="poisson")
  {
    theta_0<-0; theta_1<-0

    in_b0<- (1:length(b0))[b0 != 0]
    in_b1<- (1:length(b1))[b1 != 0]

    if( length(in_b0) !=0)
    {
      bt0<- b0[in_b0];  bt0<- as.matrix(bt0, ncol=1)
      theta_0<- A[,in_b0] %*% bt0
    }

    if( length(in_b1) !=0)
    {
      bt1<- b1[in_b1];  bt1<- as.matrix(bt1, ncol=1)
      theta_1<- A[,in_b1] %*% bt1
    }

    theta_t<- pmax(theta_0, theta_1)

    bb<- b1-b0
    in_b<- (1:length(bb))[bb != 0]

    bt<-0

    if( length(in_b) != 0)
    { bt<- bb[in_b];  bt<- as.matrix(bt, ncol=1) }

    tm1 <-  sum(bt^2) / uh
    tm2 <-  sum( exp(theta_t) * (theta_0 - theta_1)^2 )

    Bt <- exp(tm1 - tm2)
  }

  return(Bt)
}

#####Preprocessing Data
PrePro<-function(X,standardize){
  if(is.na(X)==T){
    na_index<-(1:dim(X)[2])[apply(X,2,is.na)]
    X<-X[,-na_index]
  }
  if(standardize == T){
    X<-Standardize(X)
  }
  return(as.matrix(X))
}

####likelihood function with known sigma
lh<-function(Y, X, beta, family=c("gaussian","binomial","poisson"))
{
  X<-as.matrix(X)
  n<-dim(X)[1]
  p<-dim(X)[2]

  ind_0<-(1:p)[beta != 0]
  bc_0<- as.matrix(beta[ind_0],ncol=1)
  Xs_0<-X[, ind_0]
  R_0<- Xs_0 %*% bc_0
  sig <- sum((Y-R_0)^2)/n

  return(ll<-switch(family,
            "gaussian"=-n/2*log(sig),
            "binomial"=sum( (Y*R_0 -  log(1 + exp(R_0))/sig^2 )),
            "poisson" =sum( dpois(Y, exp(R_0), log = TRUE))
            )
  )
}
####hard thresholding function


Hard<-function(t=c(0,0), k=-1, lam=1)
{

  y<-t
  t<-abs(t)

  if(k > 0)
  { lam<-sort(t, decreasing=TRUE)[k] }

  y[t<lam]<-0
  return(y)
}

#######

#####Group hard thersholding function
CI2DI<-function(I,CI_index){
  DI_index<-c();j=1;L=0
  for(i in 1:length(I$CM)){
    if(i %in% CI_index && i %in% I$CI){DI_index<-append(DI_index,i+L+seq(I$nlevel[j])-1)}
    if(i %in% CI_index && !(i %in% I$CI)){DI_index<-append(DI_index,i+L)}
    if(i %in% I$CI){L<-L+I$nlevel[j]-1;j<-j+1}
  }
  DI_index+1#for intercept
}
DI2CI<-function(I,DI_index){
  CI_index<-c();j=0;L=0
  for(i in 2:dim(I$IM)[2]){
    if(i %in% DI_index && !(i %in% unlist(I$DFI))){CI_index<-append(CI_index,i-L)}
    if(i %in% DI_index && i %in% I$DI){CI_index<-append(CI_index,i-L)}
    if(i %in% DI_index && i %in% unlist(I$DFI) && !(i %in% I$DI)){CI_index<-append(CI_index,I$DI[j]-(L-I$nlevel[j]+1))}
    if(i %in% I$DI){L<-L+I$nlevel[j+1]-1;j<-j+1}
  }
  CI_index<-unique(CI_index)
  CI_index<-CI_index-1#intercept
  CI_index
}

GroupHard<-function(Beta_t,I,k,Screening_index,penalize_mod){
  x<-Group_Beta(Beta_t,I,penalize_mod)
  CSI<-(1:length(I$CM))[Hard(t=x[Screening_index],k)==0]
  DSI<-CI2DI(I,CSI)
  Beta_t[DSI]<-0
  Beta_t
}

Group_Beta<-function(Beta_t,I,penalize_mod){
  num<-as.vector(Beta_t)
  num<-num[-unlist(I$DFI)]
  num<-num[-1]
  x<-c()
  ci<-I$CI
  for(i in 1:length(I$CM))
  {
    if(i %in% ci){
      j<-match(i,ci)
      if(penalize_mod){x<-append(x,sqrt(crossprod(Beta_t[I$DFI[[j]]]))/sqrt(I$nlevel[j]))}
      else{x<-append(x,sqrt(crossprod(Beta_t[I$DFI[[j]]])))}
    }else{
      x<-append(x,num[1])
      num<-num[-1]
    }
  }
  names(x)<-names(I$CM)
  x
}
sub_off<-function(all,sub){
  S<-all[!(all %in% sub)]
  return(S)
}
Gen_DesignMatrix<-function(correalation=c('ID','AR','MA','CS'),N=200,P=1000,index=0,rho=0.5,intercept= FALSE)
{
  this.call<-match.call()
  correalation<-match.arg(correalation)
  if (N < 0 || P < 0 || N%%1!=0 || P%%1!=0)
    stop("N,P should be postive integers")

  if(correalation=='ID'){
    X<-matrix(rnorm(N*P), nrow=N, ncol=P)
  }
  else if(correalation=='MA'){
    Z<-matrix(rnorm(N*(P+2)), nrow=N, ncol=P+2)
    X<-matrix(0, nrow=N, ncol=P)

    for(j in 1:P)
    { X[,j]<- Z[,j+2] +  Z[,j+1] + Z[,j] }

    X<-X/sqrt(3)
  }else if(correalation=='AR'){
      mat <- diag(P)
      CC<- rho^abs(row(mat)-col(mat))
      X<- rmnorm(n=N, mean=rep(0,P), varcov=CC)
  }else{
    CC<-matrix(0.3, ncol=P, nrow=P)
    CC[index, index]<- 0.15
    diag(x=CC)<-1
    X<- rmnorm(n=N, mean=rep(0,P), varcov=CC)
  }
  if(intercept){
    X<-cbind(matrix(1,nrow=dim(X)[1],ncol=1),X)
  }

  return(X)
}
Numeric2Cate<-function(data){
  Q<-quantile(data,c(0.25,0.5,0.75,1))
  for(i in 1:length(data)){
    if(as.numeric(data[i]) <Q[1]){
      data[i]<-"A"
    }else if(as.numeric(data[i])<=Q[2]){
      data[i]<-"T"
    }else if(as.numeric(data[i])<=Q[3]){
      data[i]<-"C"
    }else{data[i]<-"G"}
  }
  return(data)
}
