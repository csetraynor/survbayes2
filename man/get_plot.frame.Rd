% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_prepare.R
\name{get_plot.frame}
\alias{get_plot.frame}
\title{Prepare plot frame}
\usage{
get_plot.frame(post, strata, obs = ic2surv)
}
\arguments{
\item{post}{posterior draws in matrix format \cr}

\item{strata}{strata}

\item{obs}{original dataframe defaults to ic2surv}
}
\value{
d a plot dataset
}
\description{
This function prepares plot frame
}
