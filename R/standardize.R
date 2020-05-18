Standardize <- function(X) {
    center = colMeans(X)
    X.c = sweep(X, 2, center)
    unit.var = sqrt(apply(X.c, 2, crossprod))
    if(sum(unit.var==0)!=0){
    zero_var = (1:dim(X)[2])[unit.var==0]
    X.c<-X.c[,-zero_var]
    unit.var<-unit.var[-zero_var]
    }
    val = sweep(X.c, 2, unit.var, "/")
    return(as.matrix(val))
}
