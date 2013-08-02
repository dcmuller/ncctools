\name{ncc_sample}
\alias{ncc_sample}
\title{Generate a nested case-control study}
\usage{
  ncc_sample(entry = 0, exit, fail, origin = 0,
    controls = 1, match = list(), include = list(),
    data = NULL, silent = FALSE)
}
\arguments{
  \item{entry}{time of entry to follow-up}

  \item{exit}{time of exit from follow-up}

  \item{fail}{indicator of status on exit from follow-up
  (censored=0, fail=1)}

  \item{origin}{the origin of the analysis time-scale. For
  instance, date of birth, if age is the desired
  time-scale.}

  \item{controls}{the number of controls to sample for each
  failure}

  \item{match}{a list of categorical variables for matching
  cases and controls}

  \item{include}{a list of variables from the cohort
  dataset to be carried accross the the nested case-control
  dataset}

  \item{silent}{if \code{FALSE}, echoes a fullstop to the
  screen as each risk set is generated, and as each
  case-control set is sampled. If TRUE, suppresses output
  to the terminal.}
}
\value{
  A \code{data.frame} comprising: \item{ncc_id}{a unique
  individual identifier} \item{ncc_set}{a case-control set
  identifier} \item{ncc_fail}{case identifier (0=control,
  1=case)} \item{ncc_elig_co}{a count of the number of
  controls eligible for selection in the set}
  \item{ncc_time}{failure time of the case in the set}
  \item{ncc_pr}{the probability of being selected in the
  nested case-control sample}

  followed by the variables specified in the \code{match}
  and \code{include} lists.
}
\description{
  \code{ncc_sample} generates a nested case-control study
  dataset from a cohort study dataset. Given time of entry,
  time of exit, and exit status, risk sets are computed at
  each failure time. Controls are randomly sampled from
  these risk sets. If matching variables are specified,
  \code{ncc_sample} creates a matched or stratified nested
  case control study, in which risk sets are computed
  separately within matching strata. \code{ncc_sample} is
  similar to \code{ccwc} from the 'Epi' package, but
  differs in several small but important ways (see
  details).
}
\details{
  Given follow-up information from a cohort study,
  \code{ncc_sample} generates risk sets at each observed
  failure time, and randomly samples controls from these
  risk sets without replacement. Functionality is much the
  same as \code{ccwc} from the 'Epi' package, with two
  minor differences. Firstly, \code{ncc_sample} also
  computes and returns the total number of eligible
  controls for each risk set, as well as the probability of
  selection in to the sample for every selected case and
  control. The latter is calculated according to the
  formula given by Samuelsen (1997). Secondly,
  \code{ncc_sample} splits tied failure times at random,
  whereas \code{ccwc} preserves the ties and returns a
  multi-case case-control set.

  Control selection is performed using the
  \code{\link{sample}} function. It is therefore important
  to set the seed for the pseudo random number generator
  (\code{\link{set.seed}}) to ensure reproducibility.
}
\author{
  David C Muller
}
\references{
  Samuelsen S. O. (1997). A psudolikelihood approach to
  analysis of nested case-control studies. Biometrika,
  84(2), 379-394. doi:10.1093/biomet/84.2.379
}

