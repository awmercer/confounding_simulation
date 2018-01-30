library(bindata)
library(doRNG)

# Generate correlation matrix
cormat = function(r_xy, r_uy, r_xs, r_us, r_xu) {
  
  # Calculate correlation between S and Y based on
  # path rules
  r_ys = r_us*r_uy + r_xs*r_xy + r_us*r_xu*r_xy + r_xs*r_xu*r_uy
  
  cm = diag(4)
  rownames(cm) = c("x", "u", "y", "s")
  colnames(cm) = c("x", "u", "y", "s")
  
  cm["x", "y"] = r_xy
  cm["u", "y"] = r_uy
  cm["x", "s"] = r_xs
  cm["u", "s"] = r_us
  cm["x", "u"] = r_xu
  cm["y", "s"] = r_ys
  
  cm[lower.tri(cm)] = t(cm)[lower.tri(cm)]
  return(cm)
}

# Error handling versions of functions
quiet_rmvbin = quietly(rmvbin)
quiet_commonprob2sigma = quietly(commonprob2sigma)
quiet_check.commonprob = quietly(check.commonprob)


# Convert binary correlations to covariance matrix
# and check for validity
find_valid_sigmas = function(bincorr, margprob) {
  cp = bincorr2commonprob(margprob, bincorr)
  cp_check = quiet_check.commonprob(cp)
  if (!cp_check$result) {
    return("Invalid probabilities")
  }
  sigma = quiet_commonprob2sigma(cp)
  if (sigma$output == "") {
    rownames(sigma$result) = rownames(bincorr)
    colnames(sigma$result) = colnames(bincorr)
    return(sigma$result)
  } else {
    return("Not positive definite")
  }
}

# Create a simulated population and simulate survey samples
do_sim = function(sigma, margprob, n_pop, n_samp, n_sims, mc_cores = 1) {
  
  # Create simulated population
  pop = rmvbin(n_pop, margprob = margprob, sigma = sigma, colnames = colnames(sigma)) %>%
    as_tibble() %>%
    group_by(x, u) %>%
    mutate(pr_y_xu = mean(y),
           pr_s_xu = mean(s), 
           pi.true = pr_s_xu * n_samp/n_pop) %>%
    ungroup()
  
  popsum = pop %>%
    group_by(x, u) %>%
    summarise(n = n(),
              y_bar.true.pop = mean(y),
              pi.true = mean(pi.true) )  %>%
    ungroup()
  
  # Generate samples first so that random seed will be reproducible even if
  # simulatsions are run in parallel
  samples = rerun(n_sims, 
                  sample(seq_len(n_pop), size = n_samp, replace = FALSE, prob = pop$pr_s_xu)
  )
  
  registerDoParallel(cores = mc_cores)
  
  sim_results = foreach(i=icount(n_sims)) %dorng% {
    samp_ids = sample(seq_len(n_pop), size = n_samp, replace = FALSE, prob = pop$pr_s_xu)
    timer = timefactory::timefactory()
    
    # Identify sampled cases
    pop$samp = FALSE
    pop$samp[samp_ids] = TRUE
    
    samp_df = slice(pop, samp_ids)
    
    # Outcome regression model
    or_fit = lm(y ~ x, data = samp_df)
    
    popsamp = pop %>% 
      group_by(x, u, samp) %>% 
      summarise(n=n())
    
    # Inclusion propensity models
    pr_fit = lm(samp ~ x, data = popsamp, weights = n)

    # Estimates based on outcome regression model
    # and total population
    or_estimates = popsum %>%
      mutate(y_bar.or.pop = predict(or_fit, newdata = .),
             delta.or.pop = y_bar.or.pop - y_bar.true.pop) 
    
    # Calculate weights and other case level values for for sampled cases
    samp_df = samp_df %>%
      mutate(y_hat = predict(or_fit, newdata = samp_df), 
             pi_hat = predict(pr_fit, newdata = samp_df, type = "response"),
             wt_conf = 1/pi_hat,
             delta_y = y_hat - pr_y_xu,
             or_resid = y - y_hat)
    
    # Generate propensity weighted estimates and other values based
    # on sampled cases
    pr_estimates = samp_df %>%
      group_by(x, u) %>%
      summarise(y_bar.unwt.samp = mean(y),
                y_bar.pr.samp = weighted.mean(y, wt_conf),
                wtd_or_resid.samp = weighted.mean(or_resid, wt_conf),
                pi_hat.samp = mean(pi_hat),
                n.samp = n(),
                sum_wt.samp = sum(wt_conf))
    
    # Merge and return all estimated values
    res = left_join(or_estimates, pr_estimates, by=c("x", "u"))
    
    #cat(glue("{i}/{n_sims}: {round(timer(), 1)}\n\n"))
    
    res
    
  } %>%
    bind_rows(.id = "sim_id")
}


