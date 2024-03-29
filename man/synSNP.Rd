% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synSNP.R
\docType{data}
\name{synSNP}
\alias{synSNP}
\title{Synthetic genetic association study data set}
\format{
An object of class \code{'data.frame'} with 800 rows and 10,032 columns.
}
\usage{
data(synSNP)
}
\description{
This simulated data set consists of 10,031 genetic variants (SNPs) and a 
 continuous response variable measured on 800 individuals. The 
 genotypes were sampled from genotypic distributions derived from the  
 1000 Genomes project using the \R package \pkg{sim1000G}. The genotype 
 is coded as 0, 1, or 2 by counting the number of minor alleles (the allele 
 that is less common in the sample). The continuous response variable was 
 simulated from a normal distribution with mean that depends additively 
 on the causal SNPs.
}
\examples{
\donttest{
data(synSNP)
Y_SNP <- synSNP[,1]
X_SNP <- synSNP[,-1]
fit <- SMLE(Y = Y_SNP, X = X_SNP, k = 40)
summary(fit)
plot(fit)
}
}
\references{
The 1000 Genomes Project Consortium (2015). Global reference for human genetic variation,
  \emph{Nature}, \bold{526}(7571), 68-74.s
}
\keyword{datasets}
