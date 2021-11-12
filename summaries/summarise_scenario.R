source("common.R")

dat_dir <- "~/out_files/optimum_sims/"
fig_dir <- "figs"
tab_dir <- "tabs"
save_plots <- TRUE
scenarios <- 1:3

for(z in scenarios) {

  dat <- readRDS(paste0(dat_dir, "scenario_", z, "_dat.rds"))
  dec <- readRDS(paste0(dat_dir, "scenario_", z, "_dec.rds"))

  p_sup <- plot_superiority(dec)
  p_ess <- plot_expected_sample_size(dec)
  p_suc <- plot_expected_success(dec)
  p_fut <- plot_futility(dec)
  t_sum <- table_summary(dec)

  if (save_plots) {

    if (!dir.exists(fig_dir)) dir.create(fig_dir)
    if (!dir.exists(tab_dir)) dir.create(tab_dir)

    ggsave(paste0(fig_dir, "/scenario_", z, "_fut.png"), p_fut, width = 4.8, height = 6, dpi = 300)
    ggsave(paste0(fig_dir, "/scenario_", z, "_suc.png"), p_suc, width = 4.8, height = 6, dpi = 300)
    ggsave(paste0(fig_dir, "/scenario_", z, "_ess.png"), p_ess, width = 4.8, height = 6, dpi = 300)
    ggsave(paste0(fig_dir, "/scenario_", z, "_sup.png"), p_sup, width = 4.8, height = 6, dpi = 300)
    writeLines(t_sum, paste0(tab_dir, "/scenario_", z, "_summary.tex"))

  }

}


