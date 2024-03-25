library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(kableExtra)

theme_set(theme_bw(base_size = 9) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "grey90")))

#' @title add_facet_labs
#'
#' @param p A ggplot2
#' @param labelT Label for column facets (top)
#' @param labelR Label for row facets (right)
#' @return A ggplot2 object with facet labels added
add_facet_labs <- function(p, labelT = "", labelR = "") {
  g <- ggplotGrob(p)
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(g$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(g$layout, grepl("strip-t", name), select = t:r)
  # Add a new column to the right of current right strips,
  # and a new row on top of current top strips
  if(nrow(posR) > 0)
    width <- g$widths[max(posR$r)]    # width of current right strips
  if(nrow(posT) > 0)
    height <- g$heights[min(posT$t)]  # height of current top strips
  if(nrow(posR) > 0)
    g <- gtable_add_cols(g, width, max(posR$r))
  if(nrow(posT) > 0)
    g <- gtable_add_rows(g, height, min(posT$t)-1)

  # Construct the new strip grobs
  if(nrow(posR) > 0) {
    stripR <- gTree(name = "Strip_right", children = gList(
      rectGrob(gp = gpar(col = "black", fill = "grey90")),
      textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, fontface = 'bold', col = "grey10"))))
  }
  if(nrow(posT) > 0) {
    stripT <- gTree(name = "Strip_top", children = gList(
      rectGrob(gp = gpar(col = "black", fill = "grey90")),
      textGrob(labelT, gp = gpar(fontsize = 8.8, fontface = 'bold', col = "grey10"))))
  }

  # Position the grobs in the gtable
  if(nrow(posR) > 0) {
    g <- gtable_add_grob(g, stripR, t = min(posR$t)+1,
                         l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  }
  if(nrow(posT) > 0) {
    g <- gtable_add_grob(g, stripT, t = min(posT$t),
                         l = min(posT$l), r = max(posT$r), name = "strip-top")
  }

  # Add small gaps between strips
  if(nrow(posR) > 0)
    g <- gtable_add_cols(g, unit(1/5, "line"), max(posR$r))
  if(nrow(posT) > 0)
    g <- gtable_add_rows(g, unit(1/5, "line"), min(posT$t))
  return(g)
}


#' @title plot_superiority
#'
#' Create superiority plot from simulations
#'
#' @param decdat The decision results from simulated trials.
#' @return A ggplot2 summarising probability of superiority.
plot_superiority <- function(decdat) {
  pdat <- decdat[,
                 .(t1err = mean(fin == "success")),
                 by = .(p2tru, sup_k, fut_k, suc_k)]
  p <- ggplot(pdat,
              aes(x = sup_k, y = t1err, shape = suc_k, colour = suc_k)) +
    facet_grid(p2tru ~ fut_k, scales = "free_y") +
    # geom_point(size = 1, fill = "white") +
    geom_line() +
    scale_x_continuous(breaks = c(0.95, 0.96)) +
    scale_shape_manual(values = 0:5) +
    scale_color_viridis_d(begin = 0, end = 0.9, option = "B") +
    theme(panel.grid.minor = element_blank()) +
    labs(y = "Probability declare superiority",
         x = bquote("Superiority bound"~"("*q*")"),
         linetype = bquote("Expect success bound"~"("*bar(c)[k]*")"),
         shape = bquote("Expect success bound"~"("*bar(c)[k]*")"),
         colour = bquote("Expect success bound"~"("*bar(c)[k]*")")) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title.position="top", title.hjust = 0.5, nrow = 1),
           linetype = guide_legend(title.position="top", title.hjust = 0.5, nrow = 1),
           shape = guide_legend(title.position="top", title.hjust = 0.5, nrow = 1))
  g <- add_facet_labs(p,
                      bquote(bold(Futility)~bold(bound)~"("*underline(c)[k]*")"),
                      bquote(bold(True)~bold(probability)~"("*theta[w]^"*"*")"))
  return (grid.arrange(g))
}


#' @title plot_expected_sample_size
#'
#' Create expected sample size plot from simulations
#'
#' @param decdat The decision results from simulated trials.
#' @return A ggplot2 summarising probability of superiority.
plot_expected_sample_size <- function(decdat) {
  pdat <- decdat[,
                 .(t1err = mean(enro)),
                 by = .(p2tru, sup_k, fut_k, suc_k)]
  p <- ggplot(pdat,
              aes(x = sup_k, y = t1err, colour = suc_k, shape = suc_k)) +
    facet_grid(p2tru ~ fut_k, scales = "free_y") +
    # geom_point(size = 1) +
    geom_line() +
    scale_x_continuous(breaks = c(0.95, 0.96)) +
    scale_color_viridis_d(begin = 0, end = 0.9, option = "B") +
    theme(panel.grid.minor = element_blank()) +
    labs(y = "Expected sample size",
         x = bquote("Superiority bound"~"("*q*")"),
         colour = bquote("Expect success bound"~"("*bar(c)[k]*")"),
         shape = bquote("Expect success bound"~"("*bar(c)[k]*")")) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title.position="top", title.hjust = 0.5, nrow = 1))
  g <- add_facet_labs(p,
                      bquote(bold(Futility)~bold(bound)~"("*underline(c)[k]*")"),
                      bquote(bold(True)~bold(probability)~"("*theta[w]^"*"*")"))
  return (grid.arrange(g))
}


#' @title plot_futility
#'
#' Create futility plot from simulations
#'
#' @param decdat The decision results from simulated trials.
#' @return A ggplot2 summarising probability of superiority.
plot_futility <- function(decdat) {
  dendat <- decdat[enro < 3000 & fut_k > 0 & suc_k == 0.95,
                   .(denN = .N), keyby = .(fut_k, sup_k, p2tru)]
  pdat <- decdat[res == "futile" & suc_k == 0.95][
    enro < 3000, .N, keyby = .(stage, fut_k, sup_k, p2tru)][
      dendat, on = .(fut_k, sup_k, p2tru)][, .(P = N / denN),
                                           keyby = .(stage, fut_k, sup_k, p2tru)]
  p <- ggplot(pdat,
              aes(stage, P, colour = factor(fut_k))) +
    facet_grid(p2tru ~ sup_k, scale = "free_y") +
    geom_line() +
    scale_color_viridis_d(begin = 0, end = 0.9, option = "B") +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    ylim(0, NA) +
    labs(x = "Interim analysis (k)", y = "Marginal stopping probability",
         colour = bquote(Futility~bound~"("*underline(c)[k]*")")) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title.position="top", title.hjust = 0.5))
  g <- add_facet_labs(p,
                      bquote(bold(Superiority)~bold(bound)~("q")),
                      bquote(bold(True)~bold(probability)~"("*theta[w]^"*"*")"))
  return (grid.arrange(g))
}


#' @title plot_expected_success
#'
#' Create expected success plot from simulations
#'
#' @param decdat The decision results from simulated trials.
#' @return A ggplot2 summarising probability of superiority.
plot_expected_success <- function(decdat) {
  dendat <- decdat[enro < 3000 & fut_k == 0.05 & suc_k < 1,
                   .(denN = .N), keyby = .(suc_k, sup_k, p2tru)]
  pdat <- decdat[res == "expect success" & fut_k == 0.05][
    enro < 3000, .N, keyby = .(stage, suc_k, sup_k, p2tru)][
      dendat, on = .(suc_k, sup_k, p2tru)][, .(P = N / denN),
                                           keyby = .(stage, suc_k, sup_k, p2tru)]
  p <- ggplot(pdat,
              aes(stage, P, colour = factor(suc_k))) +
    facet_grid(p2tru ~ sup_k, scales = "free_y") +
    geom_line() +
    scale_color_viridis_d(begin = 0, end = 0.9, option = "B") +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11)) +
    ylim(0, NA) +
    labs(x = "Interim analysis (k)", y = "Marginal stopping probability",
         colour = bquote(Expect~Success~bound~"("*bar(c)[k]*")")) +
    theme(legend.position = "top") +
    guides(colour = guide_legend(title.position="top", title.hjust = 0.5))

  g <- add_facet_labs(p,
                      bquote(bold(Superiority)~bold(bound)~("q")),
                      bquote(bold(True)~bold(probability)~"("*theta[w]^"*"*")"))
  return (grid.arrange(g))
}


#' @title table_summary
#'
#' Create summary table from trial simulation
#'
#' @param decdat The decision results from simulated trials.
#' @return A summary table using kableExtra
table_summary <- function(decdat, nmax = 2000, ...) {
  tdat <- decdat[fut_k == 0.05 & suc_k == 0.95 & sup_k == 0.95,
              .("Decide superior" = mean(fin == "success"),
                "Stop early superior" = mean(fin == "success" & resp < nmax),
                "No stop superior" = mean(fin == "success" & resp == nmax),
                "Stop futile" = mean(res == "futile"),
                "Stop expect success" = mean(res == "expect success"),
                "Superior following futile" = sum(fin == "success" & res == "futile") / sum(res == "futile"),
                "Superior following expect success" =
                  sum(fin == "success" & res == "expect success") / sum(res == "expect success"),
                "Stop early" = mean(resp < nmax),
                "Expected sample size" = mean(enro)), keyby = .("$\\theta_a^\\star$" = p1tru, "$\\theta_w^\\star$" = p2tru)
  ]
  kable(tdat, "latex", digits = 2, linesep = " ", escape = FALSE, booktabs = TRUE) %>%
    kable_styling()
}
