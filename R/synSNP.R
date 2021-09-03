#' Synthetic genetic association study dataset
#'
#'  The simulated dataset consists of 10,031 genetic variants (SNPs) and a 
#'  continuous response variable measured on 800 individuals (rows). The 
#'  genotypes were sampled from genotypic distributions derived from the  
#'  1000 Genomes project, Phase 1 using the R package 
#'  sim1000G. The genotype data is coded as 0, 1, or 2 by counting the 
#'  number of minor alleles (the allele that is less common in the sample). 
#'  The continuous response variable was simulated from a normal distribution 
#'  with mean that depends additively on the causal SNPs. 
#' @references   
#' The 1000 Genomes Project Consortium (2015). Global reference for human genetic variation,
#'   \emph{Nature}, \bold{526}(7571), 68-74.
#' @docType data
#'
#' @usage data(synSNP)
#' @format An object of class \code{"matrix"}.
#' @examples
#' \donttest{
#' data(synSNP)
#' Y_SNP<- synSNP[,1]
#' X_SNP<- synSNP[,-1]
#' fit<-SMLE(Y_SNP,X_SNP,k=40)
#' }

"synSNP"