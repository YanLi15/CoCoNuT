##--------------------------------------------------------------------------------------------
## Simulate data from normal distribution for replicability analysis: one auxiliary covariate
##--------------------------------------------------------------------------------------------
#' @title Simulate data from normal distribution for replicability analysis: one auxiliary covariate
#'
#' @param m the number of features to be simulated.
#' @param xi a numeric vector of the prior probabilities for the hidden states of two primary studies corresponding to (0, 0), (0, 1), (1, 0), (1, 1).
#' @param info.str a string representing the informative strength of the auxiliary covariate: "uninormative", "weak", "moderate", "strong".
#' @param mu1 a numeric value of the signal strength in study 1.
#' @param mu2 a numeric value of the signal strength in study 2.
#' @param mu3 a numeric value of the signal strength in the auxiliary covariate.
#'
#' @return a list with the following elements:
#' \item{pvals1}{a numeric vector of the simulated p-values for study 1.}
#' \item{pvals2}{a numeric vector of the simulated p-values for study 2.}
#' \item{x}{a numeric vector of the simulated p-values for the auxiliary study.}
#' \item{theta1}{a numeric vector of 0/1 indicating the hidden states of study 1.}
#' \item{theta2}{a numeric vector of 0/1 indicating the hidden states of study 2.}
#' \item{zeta}{a numeric vector of 0/1 indicating the hidden states of the auxiliary study.}
rep_data_gen <- function(m = 10000, xi = c(0.9, 0.025, 0.025, 0.05),
                         info.str = "moderate", mu1 = 2, mu2 = 2, mu3 = 2)
{
  h = sample(0:3, m, replace = TRUE, prob = xi)
  states1 = rep(0, m)
  states1[c(which(h==2), which(h==3))] = 1
  states2 = rep(0, m)
  states2[c(which(h==1), which(h==3))] = 1

  stat1 = rnorm(m, states1*mu1, 1)
  stat2 = rnorm(m, states2*mu2, 1)

  p1 = 1 - pnorm(stat1, mean = 0, sd = 1)
  p2 = 1 - pnorm(stat2, mean = 0, sd = 1)

  truth = states1 * states2
  if(info.str == "uninformative"){
    states3 = sample(0:1, m, replace = TRUE, prob = c(sum(xi[1:2]),sum(xi[3:4])))
  }else if(info.str == "weak"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.5, 0.5))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
    # states3 = abs(truth-rand.s)
  }else if(info.str == "moderate"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.7, 0.3))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
    # states3 = abs(truth-rand.s)
  }else{
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
    # states3 = abs(truth-rand.s)
  }
  stat3 = rnorm(m, states3*mu3, 1)
  p3 = 1 - pnorm(stat3, mean = 0, sd = 1)

  return(list(pvals1 = p1, pvals2 = p2, x = p3, theta1 = states1, theta2 = states2, zeta = states3))
}

##--------------------------------------------------------------------------------------------
## Simulate data from normal distribution for replicability analysis: two auxiliary covariates
##--------------------------------------------------------------------------------------------
#' @title Simulate data from normal distribution for replicability analysis: two auxiliary covariates
#' @param m the number of features to be simulated.
#' @param xi a numeric vector of the prior probabilities for the hidden states of two primary studies corresponding to (0, 0), (0, 1), (1, 0), (1, 1).
#' @param info.str a string vector representing the informative strengths of two auxiliary covariates, with potential elements: "uninormative", "weak", "moderate", "strong".
#' @param mu1 a numeric value of the signal strength in primary study 1.
#' @param mu2 a numeric value of the signal strength in primary study 2.
#' @param mu3 a numeric value of the signal strength in the auxiliary study 1.
#' @param mu4 a numeric value of the signal strength in the auxiliary study 2.
#'
#' @return a list with the following elements:
#' \item{pvals1}{a numeric vector of the simulated p-values for primary study 1.}
#' \item{pvals2}{a numeric vector of the simulated p-values for primary study 2.}
#' \item{x1}{a numeric vector of the simulated p-values for auxiliary study 1.}
#' \item{x2}{a numeric vector of the simulated p-values for auxiliary study 2.}
#' \item{theta1}{a numeric vector of 0/1 indicating the hidden states of primary study 1.}
#' \item{theta2}{a numeric vector of 0/1 indicating the hidden states of primary study 2.}
#' \item{zeta1}{a numeric vector of 0/1 indicating the hidden states of auxiliary study 1.}
#' \item{zeta2}{a numeric vector of 0/1 indicating the hidden states of auxiliary study 2.}
rep_data_gen_multi <- function(m = 10000, xi = c(0.9, 0.025, 0.025, 0.05),
                               info.str = c("moderate", "moderate"), mu1 = 2,
                               mu2 = 2, mu3 = 2, mu4 = 2)
{
  h = sample(0:3, m, replace = TRUE, prob = xi)
  states1 = rep(0, m)
  states1[c(which(h==2), which(h==3))] = 1
  states2 = rep(0, m)
  states2[c(which(h==1), which(h==3))] = 1

  stat1 = rnorm(m, states1*mu1, 1)
  stat2 = rnorm(m, states2*mu2, 1)

  p1 = 1 - pnorm(stat1, mean = 0, sd = 1)
  p2 = 1 - pnorm(stat2, mean = 0, sd = 1)

  truth = states1 * states2

  if(info.str[1] == "uninformative"){
    states3 = sample(0:1, m, replace = TRUE, prob = c(sum(xi[1:2]),sum(xi[3:4])))
  }else if(info.str[1] == "weak"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.5, 0.5))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
    # states3 = abs(truth-rand.s)
  }else if(info.str[1] == "moderate"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.7, 0.3))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
    # states3 = abs(truth-rand.s)
  }else{
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
    states3 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states3[j] = truth[j]
      else
        states3[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
    # states3 = abs(truth-rand.s)
  }
  stat3 = rnorm(m, states3*mu3, 1)
  p3 = 1 - pnorm(stat3, mean = 0, sd = 1)

  if(info.str[2] == "uninformative"){
    states4 = sample(0:1, m, replace = TRUE, prob = c(sum(xi[1:2]),sum(xi[3:4])))
  }else if(info.str[2] == "weak"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.5, 0.5))
    states4 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states4[j] = truth[j]
      else
        states4[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
  }else if(info.str[2] == "moderate"){
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.7, 0.3))
    states4 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states4[j] = truth[j]
      else
        states4[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
  }else{
    rand.s = sample(0:1, m, replace = TRUE, prob = c(0.9, 0.1))
    states4 = rep(0, m)
    for(j in 1:m){
      if(rand.s[j]==0)
        states4[j] = truth[j]
      else
        states4[j] = sample(0:1, 1, prob = c(sum(xi[1:2]),sum(xi[3:4])))
    }
  }
  stat4 = rnorm(m, states4*mu4, 1)
  p4 = 1 - pnorm(stat4, mean = 0, sd = 1)

  return(list(pvals1 = p1, pvals2 = p2, x1 = p3, x2 = p4,
              theta1 = states1, theta2 = states2, zeta1 = states3, zeta2 = states4))
}


