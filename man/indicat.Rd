% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/qfunction.R
\name{indicat}
\alias{indicat}
\title{find out indicator for a2 for each individual}
\usage{
indicat(outcome, tx1, tx2, cov1, cov2)
}
\arguments{
\item{outcome}{is CD4 # at the end}

\item{tx1,tx2}{are treatment choices at each decision point}

\item{cov1,cov2}{are two covariate of each individual}
}
\value{
indicator of a2, if value>0, a2=1,taking a2 treatment; if value<0, a2=0,not taking a2
}
\description{
It is the best treatment choice at time 2, which maximizes Q2 function and get V2 function.
}

