% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/qfunction.R
\name{indicat.b}
\alias{indicat.b}
\title{find out indicator for a1 for each individual}
\usage{
indicat.b(outcome, tx1, tx2, cov1, cov2)
}
\arguments{
\item{outcome}{is CD4 # at the end}

\item{tx1,tx2}{are treatment choices at each decision point}

\item{cov1,cov2}{are two covariate of each individual}
}
\value{
maximum value of outcome
}
\description{
It is the best treatment choice at time 1, which maximizes Q1 function and get V1 function.
}

