% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_coxfit.R
\name{mod_coxfit}
\alias{mod_coxfit}
\title{Model fitting}
\usage{
mod_coxfit(x, surv_form, iter = 20, inits = NA_character_, ...)
}
\arguments{
\item{x}{data}

\item{mod}{Coxph model object fitted with coxph (survival).}
}
\value{
mod
}
\description{
Fits coxph model.
}
\references{
Terry M. Therneau and Patricia M. Grambsch (2000).
  Modeling Survival Data: Extending the Cox Model.
  Springer, New York. ISBN 0-387-98784-3.
}
\seealso{
[coxph]
}
\author{
Carlos S Traynor
}
\keyword{coxph}
