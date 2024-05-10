############################################################################################
# Package: coconut
# Version: 0.1.0
# Data: 2024-05-09
############################################################################################
#' @title Covariate-assisted composite null hypothesis testing providing FDR control
#'
#' @param pa a numeric vector of p-values from primary study 1.
#' @param pb a numeric vector of p-values from primary study 2.
#' @param x a numeric vector of the integrated auxiliary covariate.
#'
#' @return a list with the following elements:
#' \item{call}{the function call made.}
#' \item{Lfdr}{a numeric vector of the estimated values of joint local FDR.}
#' \item{radj}{a numeric vector of the adjusted p-values for FDR control based on the joint local FDR.}
#' \item{xi}{a numeric vector of the estimated probabilities for the joint hidden state of study 1, study 2 and auxiliary covariate.}
#' \item{f1}{a vector of the estimated densities for the p-value under the non-null in study 1}
#' \item{f2}{a vector of the estimated densities for the p-value under the non-null in study 1}
#' \item{f2}{a vector of the estimated densities for the p-value under the non-null in the auxiliary covariate.}
#' \item{loglik}{a numeric vector of the log likelihood at each step}
#'
#' @author Yan Li, \email{li_yan@@cust.edu.cn}
#' @references CoCoNuT: covariate-assisted composite null hypothesis testing with applications to replicability analysis of high-throughput experiments
#' @importFrom qvalue pi0est
#'
#' @examples
#' ## obtained two sequence of primary p-values and two sequences of auxiliary covariates: pa, pb, x1, x2
#' data.obj <- rep_data_gen_multi(m = 10000, xi = c(0.9, 0.025, 0.025, 0.05), info.str = c("moderate", "moderate"), mu1 = 2, mu2 = 2, mu3 = 2, mu4 = 2)
#' pvals1 = data.obj$pvals1
#' pvals2 = data.obj$pvals2
#' x1 = data.obj$x1
#' x2 = data.obj$x2
#'
#' ## replicability analysis exploiting the auxiliary covariates
#' x = cauchy(rbind(x1, x2))
#' res.coconut <- coconut(pvals1, pvals2, x)
#' coco.radj <- res.coconut$radj
#'
#' ## calculate the false discovery proportion and true positive proportion
#' states1 = data.obj$theta1
#' states2 = data.obj$theta2
#' fdr = sum(coco.radj <= 0.05 & !(states1 * states2))/ max(sum(coco.radj <= 0.05), 1)
#' power = sum(coco.radj <= 0.05 & (states1 * states2)) / sum(states1 * states2)
#'
#' @export
#'
coconut <- function(pa, pb, x){
  require(qvalue)
  pvals.cutoff = 1e-15
  pa[pa == 0] <- min(min(pa[pa != 0]), pvals.cutoff)
  pb[pb == 0] <- min(min(pb[pb != 0]), pvals.cutoff)
  x[x==0] <- min(min(pb[pb != 0]), pvals.cutoff)

  pi0_pa <- min(pi0est(pa)$pi0, 0.999)
  pi0_pb <- min(pi0est(pb)$pi0, 0.999)
  pi0_x <- min(pi0est(x)$pi0, 0.999)

  res <- em_pava(pa, pb, x, pi0_pa, pi0_pb, pi0_x)

  return(res)
}
