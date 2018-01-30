library(tidyverse)
library(glue)
library(doParallel)
library(DBI)
library(pool)

pool = dbPool(RSQLite::SQLite(), dbname = "data/simdata.sqlite")
# if (dbExistsTable(pool, "sim_results")) {
#   db_drop_table(pool, "sim_results")
# }

source("code/sim_functions.R")


corr_vals = seq(-.3, .3, .15)
marg_props = c(.5, .5, .5, .5)
var_labels = c("x", "u", "y", "s")

# Combinations of correlations to test
corr_combinations = cross_df(list(r_xy = corr_vals,
                                  r_uy = corr_vals,
                                  r_xs = corr_vals,
                                  r_us = corr_vals,
                                  r_xu = 0)) %>%
  filter(r_us * r_uy != 0) # Exclude cases where r_us or r_uy is 0


# Generate correlation matrices for each combination
corr_matrices = pmap(corr_combinations, cormat) 

# Convert correlation matrices to covariance matrices.
sigmas = corr_matrices %>% map(find_valid_sigmas, margprob = marg_props)

# Identify the ones that are invalid
corr_combinations$valid = map_lgl(sigmas, negate(is.character))

# Filter down to valid sigmas/correlations
sigmas = keep(sigmas, negate(is.character))
corr_combinations = filter(corr_combinations, valid) %>% 
  select(-valid) %>%
  mutate(corr_id = 1:n())

n_pop = 500000
n_samp = 500
n_sims = 100
margprob = marg_props
mc_cores = 13

set.seed(123)
t = timefactory::timefactory()
sim_results = imap(sigmas, function(sigma, corr_id) {
  cat(glue("Starting sim number {corr_id}/{length(sigmas)}\n\n"))
  t1 = timefactory::timefactory()
  res = do_sim(sigma = sigma, 
               margprob = marg_props, 
               n_pop = n_pop, 
               n_samp = n_samp, 
               n_sims = n_sims, 
               mc_cores = mc_cores) %>%
    mutate(corr_id = corr_id)
  dbWriteTable(pool, "sim_results", res, append=TRUE)
  cat(glue("Finished im number {corr_id}/{length(sigma)} in {round(t1(), 1)} seconds.\n\n"))
  invisible(res)
})
cat(glue("Finished {n_sims} simulations in {round(t(), 1)} seconds.\n\n"))












