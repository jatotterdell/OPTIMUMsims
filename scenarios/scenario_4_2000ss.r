# Scenario 1 - planned accrual of 16 participants per week ----

# Packages ----
library(OPTIMUMsims)
library(data.table)
library(parallel)

# Simulation parameters ----
sims <- 10000
cores <- 15
stage_n <- seq(200, 2000, 200)
p1tru <- rep(0.1, 6)
p2tru <- seq(0.1, 0.05, -0.01)
resp_delay <- function(n) runif(n, 48L, 72L)
ppos_q <- c(0.95, 0.955, 0.96, 0.965, 0.97)

bounds <- expand.grid(f = c(0, 0.025, 0.05, 0.075, 0.1, 0.2),
                      s = c(0.8, 0.9, 0.925, 0.95, 0.975, 1),
                      q = c(0.95,0.955,0.96,0.965, 0.97))
bounds$bounds_id <- seq_len(nrow(bounds))
bounds <- data.table(bounds, key = 'bounds_id')

# Generate data ----
dat <- rbindlist(lapply(seq_len(length(p2tru)), function(k) {
  rbindlist(mclapply(seq_len(sims), function(z) {
    est_trial_prob(
      agg_trial_dat(
        sim_trial_dat(
          p1tru = p1tru[k],
          p2tru = p2tru[k],
          nmax = max(stage_n),
          enro_rate = 5,
          resp_delay = resp_delay,
          simple_rand = FALSE
      ),
      stage_n = stage_n, min_rem = 100),
    ppos_q = ppos_q)
}, mc.cores = cores), idcol = "sim_id")}))

# Decide trial using specified bounds ----
split_dat <- split(dat, f = list(dat[, sim_id], dat[, p2tru]))
dec <- rbindlist(mclapply(seq_len(nrow(bounds)), function(i)  {
  rbindlist(lapply(split_dat, function(x) dec_trial(
    x,
    sup_k = bounds[i, q],
    fut_k = bounds[i, f],
    suc_k = bounds[i, s])))
}, mc.cores = cores), idcol = "bound_id")

# Save results ----
if (!dir.exists("~/out_files/optimum_sims")) dir.create("~/out_files/optimum_sims")
saveRDS(dat, "~/out_files/optimum_sims/scenario_4_dat.rds")
saveRDS(dec, "~/out_files/optimum_sims/scenario_4_dec.rds")
