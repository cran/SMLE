## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(SMLE)

set.seed(1)

Data_eg <- Gen_Data(n = 400, p = 1000, family = "gaussian",
correlation = "AR", rho=0.5,pos_truecoef = c(1,3,5,7,9),
effect_truecoef = c(0.5,1,-1,2,-2))

Data_eg


## -----------------------------------------------------------------------------
fit1 <- SMLE(Y = Data_eg$Y,X = Data_eg$X, k = 10, family = "gaussian")
summary(fit1)

## -----------------------------------------------------------------------------
coef(fit1)

## ----eval=F-------------------------------------------------------------------
# plot(fit1)

## -----------------------------------------------------------------------------
fit1_s <- smle_select(fit1, criterion = "ebic")
summary(fit1_s)

## ----fig.cap="Figure 3 - Plot after SMLE selection"---------------------------
plot(fit1_s)

## -----------------------------------------------------------------------------
set.seed(1)
Data_sim2 <- Gen_Data(n = 700, p = 2000, family = "binomial", num_ctgidx = 5, 
                      pos_ctgidx = c(1,3,5,7,9), effect_truecoef= c(2,2,2,-1,0.5),
                      pos_truecoef = c(1,3,5,8,12), level_ctgidx = c(3,3,3,4,5))
training=sample(c(rep(TRUE,500),rep(FALSE,200)))

Y=Data_sim2$Y
dat=data.frame(Y,Data_sim2$X)

traindat=subset(dat,training)
testdat=subset(dat,!training)

## -----------------------------------------------------------------------------
fit_1 <- SMLE(Y ~ . , family = "binomial", k = 10, data=traindat)
fit_1

## -----------------------------------------------------------------------------
threshold=0.5
values=1*(predict(fit_1, newdata = testdat,type="response")>threshold)
(confusionmat=table(Actual=testdat$Y,Predicted=values))
(accuracy=sum(diag(confusionmat))/sum(confusionmat))
(ppv=confusionmat[2,2]/colSums(confusionmat)[2])
(npv=confusionmat[1,1]/colSums(confusionmat)[1])

