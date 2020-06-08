#tools function for selection.
#Multicores implenments for features selection.
ebicc<-function(Y,X_s,family,tune,k_min,k_max,

                n,pp,gamma,para,num_cores){
  #fast_smle_by_sparisty function
  #return retained feature ids for given sparisity
  fss<-function(v){

    ff<- SMLE(Y=Y, X=X_s, k=v, family=family,fast=TRUE)

    return(lh(Y=Y, X=X_s[,ff$ID_Retained],
              beta=as.vector(glm(Y ~ X_s[,ff$ID_Retained]-1, family=family)$coefficients),
              family=family))
  }

  ll<-unlist(lapply(k_min:k_max,fss))

  select_crit<-switch(tune,
                      'ebic'=function(ll,v){return(-2 * ll  + v* log(n) + 2 * v * gamma * log(pp))},
                      'bic' =function(ll,v){return(-2 * ll +  v* log(n))},
                      'aic' =function(ll,v){return(-2 * ll +  v*2)}
                      )
  return(mapply(select_crit,ll,k_min:k_max))
}


ctg_ebicc<-function(Y,X_s,family,tune,codingtype,

                   k_min,k_max,n,pp,gamma,para,num_cores){

  #fast_smle_by_sparisty function



  #return retained feature ids for given sparisity

  fss<-function(v,codingtype){

    ff<- SMLE(Y=Y, X=X_s, k=v, family=family,codingtype = codingtype)

    X_sub <- X_s[,ff$ID_Retained]

    X_dummy <- suppressWarnings(dummy.data.frame(X_sub ,sep="_",codingtype = codingtype))

    Ci <- sapply(X_sub,is.factor)

    if(codingtype =="all"){

      dummy_f <- sum(sapply(X_sub[,Ci],nlevels)-1)

    }else{

      dummy_f <- sum(sapply(X_sub[,Ci],nlevels)-2)

    }

    ll <- lh(Y=Y, X=as.matrix(X_dummy, ncol= dim(X_dummy)[2]),
             beta=as.vector(glm(Y ~ as.matrix(X_dummy,ncol=dim(X_dummy)[2])-1, family=family)$coefficients),
             family=family)
    message(ll)
    return(list(d_f= dummy_f, likelihood = ll))

  }

  select<-lapply(k_min:k_max,fss,codingtype=codingtype)

  select_crit<-switch(tune,
                      'ebic'=function(select,v){return(-2 *  select[[v]][[1]]  +   select[[v]][[2]]* log(n) + 2 *  select[[v]][[2]] * gamma * log(pp))},
                      'bic' =function(select,v){return(-2 *  select[[v]][[1]]  +   select[[v]][[2]]* log(n))},
                      'aic' =function(select,v){return(-2 *  select[[v]][[1]]  +   select[[v]][[2]]*2)}
  )


  return(unlist(lapply(k_min:k_max,select_crit,select=select)))
}
