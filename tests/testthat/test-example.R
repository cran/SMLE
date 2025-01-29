library(testthat)
library(SMLE)

Data_ctg <- Gen_Data(n = 200, p = 1000, family = "gaussian", pos_ctgidx = c(1, 2, 3), 
                     level_ctgidx = c(3, 4, 5))
head(Data_ctg$X)[, 1:5]


# Section 3.2---Demo code---Joint feature screening--------------------------------

fit <- SMLE(Y = Data_ctg$Y, X = Data_ctg$X, k = 15, family = "gaussian", keyset = c(1, 4, 5), 
            categorical = TRUE, group = TRUE)
fit

# Section 3.2---Demo code---Post-screening Selection-------------------------------

fit_s <- smle_select(fit, criterion = "ebic", gamma_seq = seq(0,1,0.2), vote = TRUE)
fit_s
