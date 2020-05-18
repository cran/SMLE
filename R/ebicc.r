#tools function for selection.
#Multicores implenments for features selection.
ebicc<-function(Y,X,family,tune,k_min,k_max,n,pp,gamma,para,num_cores){
  #fast_smle_by_sparisty function
  #return retained feature ids for given sparisity
  fss<-function(v){
    ff<- SMLE(Y=Y, X=X, k=v, family=family,fast=TRUE)
    return(lh(Y=Y, X=as.matrix(X[,ff$Retained_Feature_IDs],ncol=v),
              beta=as.vector(glm(Y ~ X[,ff$Retained_Feature_IDs]-1, family=family)$coefficients),
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


