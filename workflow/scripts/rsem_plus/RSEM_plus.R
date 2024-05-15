#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Combine counts.csv and rsem.csv files to perform
## transcript-isoform level analysis

### Load dependencies

library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(optparse)
library(MASS)


# Helper functions that I will use on multiple occasions
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) exp(x)/(1+exp(x))

# Likelihood function for mixture model
mixed_lik <- function(pnew, pold, TC, nT, n, logit_fn, p_sd = 1, p_mean = 0){
  logl <- sum(n*(log(inv_logit(logit_fn)*dbinom(TC, size = nT, prob = pnew) + (1-inv_logit(logit_fn))*dbinom(TC, nT, pold) ) ) ) + log(stats::dnorm(logit_fn, mean = p_mean, sd = p_sd))
  return(-logl)
}


# Density function for a skewed normal taken from sn package to reduce dependencies
dsn <- function (x, xi = 0, omega = 1, alpha = 0, tau = 0, dp = NULL, log = TRUE){

    if (!is.null(dp)) {
        if (!missing(alpha)) 
            stop("You cannot set both 'dp' and component parameters")
        xi <- dp[1]
        omega <- dp[2]
        alpha <- dp[3]
        tau <- if (length(dp) > 3) 
            dp[4]
        else 0
    }
    
    za <- cbind((x - xi)/omega, alpha)
    z <- za[, 1]
    alpha <- za[, 2]
    logN <- (-log(sqrt(2 * pi)) - logb(omega) - z^2/2)
    logS <- numeric(length(z))
    ok <- (abs(alpha) < Inf)
    logS[ok] <- pnorm(tau * sqrt(1 + alpha[ok]^2) + (alpha * 
        z)[ok], log.p = TRUE)
    logS[!ok] <- log(as.numeric((sign(alpha) * z)[!ok] + tau > 
        0))
    logPDF <- as.numeric(logN + logS - pnorm(tau, log.p = TRUE))
    logPDF <- replace(logPDF, abs(x) == Inf, -Inf)
    logPDF <- replace(logPDF, omega <= 0, NaN)
    out <- if (log) 
        logPDF
    else exp(logPDF)
    names(out) <- names(x)
    return(out)

}


### Process arguments

args = commandArgs(trailingOnly = TRUE)


### Read parameters

option_list <- list(
    make_option(c("-c", "--counts", type="character"),
                    default = ".",
                    help = "count.csv file path"),
    make_option(c("-r", "--rsem", type="character"),
                    default = ".",
                    help = 'rsem.csv file path'),
    make_option(c("-o", "--output", type = "character"),
                    default = ".",
                    help = 'output file path'),
    make_option(c("-s", "--sample", type = "character"),
                    default = "",
                    help = "sample name"),
    make_option(c("-e", "--echocode", type="logical"),
                    default = "FALSE",
                    help = 'print R code to stdout'),
    make_option(c("-n", "--pnew", type="double"),
                    default = -1,
                    help = 'mutation rate in new reads'),
    make_option(c("-b", "--pold", type="double"),
                    default = -1,
                    help = 'background mutation rate in old reads'),
    make_option(c("-i", "--niter", type="integer"),
                default = 14,
                help = 'Number of EM iterations'),
    make_option(c("-t", "--threshold", type="double"),
                default = 1e-3,
                help = 'change in log-likelihood threshold for convergence'),
    make_option(c("-w", "--priornew", type = "double"),
                default = 1,
                help = "Prior on the number of new reads from a transcript."),
    make_option(c("-p", "--priortot", type = "double"),
                default = 2,
                help = "Prior on the total number of reads from a transcript."),
    make_option(c("--priorpnewmean", type = "double"),
                default = -2,
                help = "Prior mean on the logit(pnew)."),
    make_option(c("--priorpnewsd", type = "double"),
                default = 0.8,
                help = "Prior sd on the logit(pnew)."),
    make_option(c("--priorpnewskew", type = "double"),
                default = -3,
                help = "Skew on logit(pnew) prior; technically alpha in sn::dsn."),
    make_option(c("--priorpoldmean", type = "double"),
                default = -5,
                help = "Prior mean on the logit(pnew)."),
    make_option(c("--priorpoldsd", type = "double"),
                default = 1.25,
                help = "Prior sd on the logit(pnew)."),
    make_option(c("--priorpoldskew", type = "double"),
                default = -10,
                help = "Skew on logit(pnew) prior; technically alpha in sn::dsn."))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.

options(echo = as.logical(opt$echocode))

### Check that pnew and pold make sense
if(opt$pnew != -1 & opt$pnew < 0){
  stop("pnew must be > 0 or equal to -1! Set to -1 if you want RSEM+ to estimate pnew for you.")
}

if(opt$pold != -1 & opt$pold < 0){
  stop("pold must be > 0 or equal to -1! Set to -1 if you want RSEM+ to estimate pold for you.")
}


### Load csvs 
counts <- fread(opt$counts)

# Print column names to debug
colnames(counts)
nrow(counts)

rsem <- fread(opt$rsem)
setkey(rsem, qname)

### Only keep relevant columns in counts.csv
    # Should eventually allow other mutation types to be analyzed
    # according to mutType argument
counts <- counts[,c("GF", "qname", "TC", "nT")]
counts <- counts[complete.cases(counts),]
setkey(counts, qname)

# Add RSEM's P(transcript) to mutational information
  # If pt is 0, don't keep it; would lead to a 0/0 situation in the EM step
cT <- rsem[counts, nomatch = NULL]
cT <- cT[pt > 0,]

### Estimate new and old mutation rate

if(opt$pnew == 0){

  rm(counts)
  rm(rsem)
  
  message("Provided pnew is 0, so assuming this is a -s4U control sample.")

  # Calculate logit(fn) with analytical Bayesian approach
  Fn_est <- cT[,.(fn_est = 0,
                  nreads = sum(pt)), by = .(GF, TF)]

  Fn_est$sample <- opt$sample

}else{

  cB <- counts[!grepl("__", GF), .N, by = .(GF, TC, nT)]
  mutrate_df <- cB[, .(n = sum(N)), by = .(TC, nT)]

  rm(counts)
  rm(rsem)

  if(opt$pnew == -1 & opt$pold == -1){

    ## USER PROVIDED NEITHER PNEW OR POLD

    # Binomial mixture likelihood
    mixture_lik <- function(param, TC, nT, n){

        logl <- sum(n*log(inv_logit(param[3])*dbinom(TC, nT, inv_logit(param[2])) +  (1-inv_logit(param[3]))*dbinom(TC, nT, inv_logit(param[1])) ) )

        totlogl <- -logl - dsn(param[2], 
                           opt$priorpnewmean,
                           opt$priorpnewsd,
                           opt$priorpnewskew,
                           log = TRUE)
                      - dsn(param[1],
                           opt$priorpoldmean,
                           opt$priorpoldsd,
                           opt$priorpoldskew,
                           log = TRUE) - dnorm(param[3], log = TRUE)

        return(totlogl)

    }

    low_ps <- c(-9, -9, -9)
    high_ps <- c(0, 0, 9)

    fit <- stats::optim(par=c(-7,-2,0), mixture_lik, TC = mutrate_df$TC, nT = mutrate_df$nT,
                                    n = mutrate_df$n, method = "L-BFGS-B", 
                                    lower = low_ps, upper = high_ps)


    pnew <- inv_logit(max(fit$par[1:2]))
    pold <- inv_logit(min(fit$par[1:2]))

    print(paste0("Estimated pnew is: ", pnew))
    print(paste0("Estimated pold is: ", pold))

  }else if(opt$pnew == -1 & opt$pold != -1){

    ## USER PROVIDED POLD BUT NOT PNEW

    # Binomial mixture likelihood
    mixture_lik <- function(param, TC, nT, n){

        logl <- sum(n*log(inv_logit(param[2])*dbinom(TC, nT, inv_logit(param[1])) +  (1-inv_logit(param[2]))*dbinom(TC, nT, opt$pold) ) )

        return(-logl - dsn(param[1], 
                           opt$priorpnewmean,
                           opt$priorpnewsd,
                           opt$priorpnewskew,
                           log = TRUE) - dnorm(param[2], log = TRUE))

    }

    low_ps <- c(logit(opt$pold), -9)
    high_ps <- c(0, 9)

    fit <- stats::optim(par=c(-2,0), mixture_lik, TC = mutrate_df$TC, nT = mutrate_df$nT,
                                    n = mutrate_df$n, method = "L-BFGS-B", 
                                    lower = low_ps, upper = high_ps)


    pnew <- inv_logit(fit$par[1])
    pold <- opt$pold

    print(paste0("Estimated pnew is: ", pnew))
    print(paste0("Provided pold is: ", pold))

  }else if(opt$pnew != -1 & opt$pold == -1){

    ## USER PROVIDED PNEW BUT NOT POLD


      # Binomial mixture likelihood
    mixture_lik <- function(param, TC, nT, n){

        logl <- sum(n*log(inv_logit(param[2])*dbinom(TC, nT, opt$pnew) +  (1-inv_logit(param[2]))*dbinom(TC, nT, inv_logit(param[1])) )) 

        return(-logl - dsn(param[1],
                           opt$priorpoldmean,
                           opt$priorpoldsd,
                           opt$priorpoldskew,
                           log = TRUE) - dnorm(param[2], log = TRUE))

    }

    low_ps <- c(-9, -9)
    high_ps <- c(logit(opt$pnew), 9)

    fit <- stats::optim(par=c(-7,0), mixture_lik, TC = mutrate_df$TC, nT = mutrate_df$nT,
                                    n = mutrate_df$n, method = "L-BFGS-B", 
                                    lower = low_ps, upper = high_ps)


    pnew <- opt$pnew
    pold <- inv_logit(fit$par[1])

    print(paste0("Provided pnew is: ", pnew))
    print(paste0("Estimated pold is: ", pold))

  }else{

    ## USER PROVIDED BOTH PNEW AND POLD
    
    pnew <- opt$pnew
    pold <- opt$pold

    print(paste0("Provided pnew is: ", pnew))
    print(paste0("Provided pold is: ", pold))


  }




  ### Get priors





  # Likelihood function for mixture model
  mixed_lik <- function(pnew, pold, TC, nT, n, logit_fn, p_sd = 1, p_mean = 0){
    logl <- sum(n*(log(inv_logit(logit_fn)*dbinom(TC, size = nT, prob = pnew) + (1-inv_logit(logit_fn))*dbinom(TC, nT, pold) ) ) ) + log(stats::dnorm(logit_fn, mean = p_mean, sd = p_sd))
    return(-logl)
  }

  # Estimate prior for GF
  Fn_prior <- cB %>% dplyr::ungroup() %>%
    dplyr::group_by(GF, TC, nT) %>%
    dplyr::summarise(n = sum(N)) %>%
    dplyr::group_by(GF) %>%
    dplyr::summarise(logit_fn_rep = stats::optim(0, mixed_lik, nT = nT, TC = TC, n = n, pnew = pnew, pold = pold, method = "L-BFGS-B", lower = -7, upper = 7)$par, nreads =sum(n), .groups = "keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prior = inv_logit(logit_fn_rep)) %>%
    dplyr::select(GF, prior)

  Fn_prior <- setDT(Fn_prior)

  # Add prior info
  setkey(Fn_prior, GF)
  setkey(cT, GF)
  cT <- cT[Fn_prior, nomatch = NULL]

  setkey(cT, GF, TF)

  ### Estimate fn with EM
  oldll <- -Inf

  # Likelihoods that can get pre-computed (new and old reads from transcript t)
  cT[, lnew := dbinom(TC, nT, pnew) * pt]
  cT[, lold := dbinom(TC, nT, pold) * pt]

  for(i in 1:opt$niter){    
    
    # Sum old + new likelihood for all possible isoforms
    cT[, den := sum(pt*prior*dbinom(TC, nT, pnew) + pt*(1-prior)*dbinom(TC, nT, pold)), by = qname]


    # E[# of new reads from transcript i] / E[# of reads from transcript i]
    cT[, fn_est := (sum((lnew * prior) / den) + opt$priornew) / (sum((lnew * prior + (1 - prior) * lold ) / den ) + opt$priortot),
          by = .(GF, TF)]
    cT[, prior := fn_est]

    # Calculate log-likelihood to assess convergence
      # I am doing this for each gene since convergence is independent for each gene; so this can be
      # seen as checking to ensure that all isoform estimates have converged.
      # Would be nice to iteratively filter out those genes that have already converged to improve
      # efficiency
    ll_vect <- unname(unlist(cT[,.(loglikelihood = log(sum(fn_est*lnew + (1-fn_est)*lold))), by = GF][, GF := NULL ]))
    ll_diff <- oldll - ll_vect
    
    if(max(abs(ll_diff)) < opt$threshold){
      
      print("Exiting before max number of iterations due to convergence")
      print(paste0("Number of iterations: ", i))

      break
      
    }else{
      
      oldll <- ll_vect

    }

  }

  # Get estimate for each isoform
  Fn_est <- cT[,.(fn_est = mean(fn_est)), by = .(GF, TF)]
  
  
  Fn_est[, sample := opt$sample]
}



write_csv(Fn_est, file = opt$output)
