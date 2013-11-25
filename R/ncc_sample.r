#' @title Generate a nested case-control study
#' 
#' @description
#' \code{ncc_sample} generates a nested case-control study dataset from a 
#' cohort study dataset. Given time of entry, time of exit, and exit status, 
#' risk sets are computed at each failure time. Controls are randomly sampled
#' from these risk sets. If matching variables are specified, \code{ncc_sample}
#' creates a matched or stratified nested case control study, in which risk
#' sets are computed separately within matching strata. \code{ncc_sample} is
#' similar to \code{ccwc} from the 'Epi' package, but differs in several small 
#' but important ways (see details).
#' 
#' @param entry time of entry to follow-up
#' @param exit time of exit from follow-up
#' @param fail indicator of status on exit from follow-up (censored=0, fail=1)
#' @param origin the origin of the analysis time-scale. For instance, date of
#' birth, if age is the desired time-scale.
#' @param controls the number of controls to sample for each failure
#' @param match a list of categorical variables for matching cases and controls
#' @param include a list of variables from the cohort dataset to be carried
#' accross the the nested case-control dataset
#' @param data a \code{\link{data.frame}} which contains the follow-up, 
#' matching, and included variables.
#' @param keep_all if \code{TRUE}, does not sample from the risk sets, but 
#' returns the probability of selection for each observation that is eligible 
#' to be selected for any case in \code{data}. Defaults to \code{FALSE}.
#' @param silent if \code{FALSE}, provides entertainment by echoing a fullstop 
#' to the screen as each risk set is generated. If TRUE, output to the
#' console is suppressed.
#' 
#' @details
#' Given follow-up information from a cohort study, \code{ncc_sample} generates
#' risk sets at each observed failure time, and randomly samples controls from
#' these risk sets without replacement. Functionality is much the same as 
#' \code{ccwc}  from the 'Epi' package, with two minor differences. Firstly, 
#' \code{ncc_sample} also computes and returns the total number of eligible 
#' controls for each risk set, as well as the probability of selection in to 
#' the sample for every selected case and control. The latter is calculated
#' according to the formula given by Samuelsen (1997). Secondly, 
#' \code{ncc_sample} splits tied failure times at random, whereas \code{ccwc}
#' preserves the ties and returns a multi-case case-control set.
#' 
#' Random sampling of controls within risk sets is performed using \R's 
#' pseudo-random number facilities. It is therefore important to set the seed 
#' (\code{\link{set.seed}}) to ensure reproducibility.
#' 
#' @return
#' A \code{data.frame} comprising:
#' \item{ncc_set}{a case-control set identifier}
#' \item{ncc_id}{a unique individual identifier}
#' \item{ncc_fail}{case identifier (0=control, 1=case)}
#' \item{ncc_elig_co}{a count of the number of controls eligible for selection
#' in the set}
#' \item{ncc_time}{failure time of the case in the set}
#' \item{ncc_pr}{the probability of being selected in the nested 
#' case-control sample}
#' 
#' followed by the variables specified in the \code{match} and \code{include} 
#' lists.
#' 
#' @author David C Muller
#' 
#' @references
#' Samuelsen S. O. (1997). A psudolikelihood approach to analysis of nested 
#' case-control studies. Biometrika, 84(2), 379-394. doi:10.1093/biomet/84.2.379
#' 
#' @import data.table
#' @export
#' 
ncc_sample <- function(entry = 0, exit, fail, origin = 0, controls = 1, 
                       match = list(), include = list(), data = NULL, 
                       keep_all = FALSE, silent = FALSE) {
  entry <- eval(substitute(entry), data)
  exit <- eval(substitute(exit), data)
  fail <- eval(substitute(fail), data)
  origin <- eval(substitute(origin), data)
  n <- length(fail)
# individual id for identifying sets of split records
  ncc_id <- 1:n
# add this id to the supplied data frame for later merging
  data <- data.table(data)
  data[, "ncc_id" := ncc_id, with=FALSE]
  if (length(exit) != n) {
    stop("All vectors must have length")
  }
  if (length(entry) == 1) {
    entry <- rep(entry, n)
  } 
  else {
    if (length(entry) != n) { 
      stop("All vectors must have same length")
    }
  }
  if (length(origin) == 1) {
    origin <- rep(origin, n)
  }
  else {
    if (length(origin) != n) 
      stop("All vectors must have same length")
  }
  t_entry <- as.numeric(entry - origin)
  t_exit <- as.numeric(exit - origin)
  
  # process matching vectors
  arg_match <- substitute(match)
  match <- lapply(arg_match, eval, data)
  names(match) <- as.character(arg_match)
  if (names(match[1]) == "list" && class(match[[1]]) == "function") {
    match[1] <- NULL
  }
  match_names <- names(match)
  
  # process 'included' vectors
  arg_inc <- substitute(include)
  include <- lapply(arg_inc, eval, data)
  names(include) <- as.character(arg_inc)
  if (names(include[1]) == "list" && class(include[[1]]) == "function") {
    include[1] <- NULL
  }
  inc_names <- names(include)
  
  grp <- rep(1, n)
  if (length(match) > 0) {
    grp <- as.numeric(factor(do.call(paste0, match)))
  }
  if (!silent) {
    cat("\nGenerating risk sets:\n")
  }
  fg <- unique(grp[fail != 0])
  gresult <- vector("list", length(fg))
  # index for groups
  gii <- 0
  set <- 0
  ties <- 0
  for (g in fg) {
    gii <- gii + 1
    ft <- unique(t_exit[(grp == g) & (fail != 0)])
    ft <- ft[order(ft)]
    tresult <- vector("list", length(ft))
    
    # index for times
    for (tii in 1:length(ft)) {
      t_end <- ft[tii]
      case <- (grp == g) & (t_exit == t_end) & (fail != 0)
      ncase <- sum(case)
      if (ncase > 1) {
        ties <- ties + 1 
      }
      # split tied failures
      sresult <- vector("list", ncase)
      case_index <- cumsum(case)*case
      for (tjj in 1:ncase) {
        set <- set + 1
        if (!silent) {
          dots <- ifelse((set %% 50 == 0), ". 50\n", ".")
          cat(dots)
        }
        failure <- case & case_index==tjj
        noncase <- (grp == g) & (t_entry <= t_end) & (t_exit >= t_end) & !failure
        n_eligible <- sum(noncase)
        ineligible <- !failure & !noncase
        df <- data.frame(ncc_id   = ncc_id, 
                         ncc_set  = set, 
                         ncc_fail = failure, 
                         nfail    = rep(sum(failure), length(failure)),
                         ncc_elig_co = rep(n_eligible, length(failure)), 
                         ncc_time     = rep(t_end, length(failure)),
                         ineligible   = ineligible)
        sresult[[tjj]] <- df[!ineligible, ]
      }
      tresult[[tii]] <- rbindlist(sresult)
    }
    gresult[[gii]] <- rbindlist(tresult)
    rm(tresult)
  }
  tsplit <- rbindlist(gresult)
  
# Calculate probability of inclusion 
# (Samuelsen, Biometrika, 1997, 84(2) p379)
  tsplit[, "pr_not" := 1 - ((nfail * controls) / ncc_elig_co)]
  tsplit[, "ncc_pr" := 1- prod(pr_not), by=ncc_id]
  tsplit[, "pr_not" := NULL]
# cases have probability of inclusion of 1
  tsplit[ncc_fail==TRUE, ncc_pr := 1]
  # if risk sets are smaller than requested number of cases, the calculated
# probability will be greater than 1. The correct value is 1.
  tsplit[ncc_pr > 1, ncc_pr := 1]
  # we no longer need the number of failures per set
  tsplit[, nfail := NULL] 
  
# sample controls from risk sets  
  if (!silent) {
    cat("\nSampling from risk sets:\n")
  }
#   browser()
  tsplit[, "r1" := runif(nrow(tsplit))]
  tsplit[, "r2" := runif(nrow(tsplit))]
  tsplit[, ncc_no_fail := !ncc_fail]
  setkey(tsplit, ncc_set, ncc_no_fail, r1, r2)
  # calculate conditions that should trigger warnings (incomplete sets, etc)
  tsplit[, c("nca", "nel") := list(nca = sum(ncc_fail),
                                   nel = sum(!ncc_fail)),
         by = ncc_set]
  tsplit[, nco := controls*nca]
  incomplete <- sum(tsplit[, sum(nel < nco), by = ncc_set][,V1]>0)
  nomatch <- sum(tsplit[, sum(nel == 0), by = ncc_set][,V1]>0)
  # return either the whole dataset (one record per person), or
  # just the ncc sample, depending on option keep_all
  if (!keep_all) {
  ncc_frame <- tsplit[, 
                 head(.SD, nca + nco), 
                 by = list(ncc_set)]
  }
  else {
    ncc_frame <- tsplit[, head(.SD, 1), by = list(ncc_id)]
  }
  
  # coerce failure flag to numeric 
  # (logical var for case status might confuse people)
  ncc_frame[, ncc_fail := as.numeric(ncc_fail)]
  # get rid of the interim variables
  ncc_frame[, c("ncc_elig_co", "r1", "r2", "ncc_no_fail", 
                "nca", "nel", "nco") := NULL]
  
  # merge ncc frame with original data
  setkey(ncc_frame, ncc_id)
  setkey(data, ncc_id)
  res <- data[ncc_frame]
  res <- res[, c(names(ncc_frame), match_names, inc_names), with=FALSE]
  res <- res[order(ncc_set,-ncc_fail), ]
  setkey(res, ncc_set)
  
  # display warnings and return the frame
  if (incomplete > 0) {
    warning(paste(incomplete, "case-control sets are incomplete"))
  }
  if (nomatch > 0) {
    warning(paste(nomatch, "cases could not be matched"))
  }
  if (ties > 0) {
    warning(paste("There", ifelse(ties == 1, "was", "were"), ties, 
                  ifelse(ties == 1, "set of", "sets of"),
                  "tied failure times"))
  }
  return(res)
}
  