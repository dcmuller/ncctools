% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncc_eligible_controls.r
\name{ncc_eligible_controls}
\alias{ncc_eligible_controls}
\title{Calculate the number of eligible controls for each case for a nested 
case control study}
\usage{
ncc_eligible_controls(entry = 0, exit, fail, origin = 0, match = list(),
  include = list(), data = NULL, silent = FALSE)
}
\arguments{
\item{entry}{time of entry to follow-up}

\item{exit}{time of exit from follow-up}

\item{fail}{indicator of status on exit from follow-up (censored=0, fail=1)}

\item{origin}{the origin of the analysis time-scale. For instance, date of
birth, if age is the desired time-scale.}

\item{match}{a list of categorical variables for matching cases and controls}

\item{include}{a list of variables to include in the returned data.frame}

\item{data}{a \code{\link{data.frame}} which contains the follow-up and
matching variables.}

\item{silent}{if \code{FALSE}, provides entertainment by echoing a fullstop 
to the screen as each risk set is generated. If TRUE, output to the
console is suppressed.}
}
\value{
A \code{data.frame} comprising:
\item{ncc_id}{a unique individual identifier}
\item{ncc_elig_co}{a count of the number of controls eligible for selection
in the set}

followed by the variables specified in the \code{match} and \code{include} 
lists.
}
\description{
\code{ncc_eligible_controls} calculates the number of eligible controls for each
case for a nested case control study. Given time of entry, time of exit, and 
exit status, 
risk sets are computed at each failure time. If matching variables are specified, 
\code{ncc_eligible_controls} creates matched or stratified risk sets, 
which are computed separately within matching strata.
}
\details{
Given follow-up information from a cohort study, \code{ncc_eligible controls}
 generates
risk sets at each observed failure time, and counts the number of eligible 
controls within each risk set. This information can be used to calculate the
probability of inclusion in the NCC study (see \code{\link{ncc_selection_prob}})
}
\author{
David C Muller
}
