% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary_two_arm.R
\name{calc_ppos}
\alias{calc_ppos}
\title{Calculate the predicted probability of success}
\usage{
calc_ppos(a, b, c, d, m1, m2, k_ppos, post_method = "exact", post_sim = 10000)
}
\arguments{
\item{a}{First parameter of first beta random variable}

\item{b}{Second parameter of first beta random variable}

\item{c}{First paramter of second beta random variable}

\item{d}{Second parameter of second beta random variable}

\item{m1}{Sample size to predict for first beta random variable}

\item{m2}{Sample size to predict for second beta random variable}

\item{k_ppos}{The posterior probability cut-point to be assessed}

\item{post_method}{The method to use for calculating posterior probabilities,
one of "exact" (numerical), "approx", "sim".}

\item{post_sim}{Number of posterior simulations if post_method = "sim".}
}
\value{
The predicted probability of success
}
\description{
Calculate the predicted probability of success
}
