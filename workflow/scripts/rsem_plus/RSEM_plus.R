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
                default = 7,
                help = 'Number of EM iterations'))

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

### Only keep relevant columns in counts.csv
    # Should eventually allow other mutation types to be analyzed
    # according to mutType argument
counts <- counts[,c("GF", "qname", "TC", "nT")]
counts <- counts[complete.cases(counts),]

# Add RSEM's P(transcript) to mutational information
  # If pt is 0, don't keep it; would lead to a 0/0 situation in the EM step
cT <- rsem[counts, on = .(qname), nomatch = NULL]
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

        logl <- sum(n*log(inv_logit(param[3])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[2])^TC)*((1 -inv_logit(param[2]))^(nT-TC)) +  (1-inv_logit(param[3]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 - inv_logit(param[1]))^(nT-TC)) ) )

        return(-logl)

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

        logl <- sum(n*log(inv_logit(param[2])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 -inv_logit(param[1]))^(nT-TC)) +  (1-inv_logit(param[2]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(opt$pold)^TC)*((1 - inv_logit(opt$pold))^(nT-TC)) ) )

        return(-logl)

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

        logl <- sum(n*log(inv_logit(param[2])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(opt$pnew)^TC)*((1 -inv_logit(opt$pnew))^(nT-TC)) +  (1-inv_logit(param[2]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 - inv_logit(param[1]))^(nT-TC)) ) )

        return(-logl)

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

  # Estimate GF for prior
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
  cT <- cT[Fn_prior, on = .(GF), nomatch = NULL]

  ### Estimate fn with EM
  for(i in 1:opt$niter){
    
    
    cT[,den := sum(pt*prior*dbinom(TC, nT, pnew) + pt*(1-prior)*dbinom(TC, nT, pold)), by = qname]
    
    Fn_est <- cT[,.(fn_est = sum(pt*prior*dbinom(TC, nT, pnew)/den)/sum((pt*prior*dbinom(TC, nT, pnew) + pt*(1-prior)*dbinom(TC, nT, pold))/den)), by = .(GF, TF)]
    
    Fn_est[, prior := fn_est]
    
    cT <- cT[,!c("prior")]
    cT <- cT[Fn_est, on = .(GF, TF), nomatch = NULL]
    
    
  }
  
  
  Fn_est[, sample := opt$sample]
}



write_csv(Fn_est, file = opt$output)
