#' @title Generate a nested case-control study
#' 
#' @description
#' \code{ncc_selection_prob} calculates the probability of selection into a 
#' nested case control study. Given a nested case control study with
#' time of entry, time of exit, exit status, and the number of potentially 
#' eligible controls for each case from the full cohort, the probability of 
#' being selected as a control for one or more cases is calculated.
#' 
#' @param entry time of entry to follow-up
#' @param exit time of exit from follow-up
#' @param fail indicator of status on exit from follow-up (censored=0, fail=1)
#' @param origin the origin of the analysis time-scale. For instance, date of
#' birth, if age is the desired time-scale.
#' @param controls the number of controls  that were sampled for each failure
#' @param match a list of categorical variables that were matched on
#' @param caseset a caseset identifier (common within matched case-control sets)
#' @param elig_controls an integer indicating the number of individuals from 
#' the full cohort who were eligible be selected as a control for each case. 
#' This should be constant within \code{caseset}.
#' @param data a \code{\link{data.frame}} which contains the follow-up, 
#' matching, caseset identifier, and number of eligible controls variables.
#' @param silent if \code{FALSE}, provides entertainment by echoing a fullstop 
#' to the screen as each risk set is generated. If TRUE, output to the
#' console is suppressed.
#' 
#' @details
#' @return
#' A \code{data.frame} comprising the original data, with a variable indicating
#' the probability of being selected in the nested case control sample (ncc_pr)
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
ncc_selection_prob <- function(entry = 0, exit, fail, origin = 0, controls = 1, 
                       match = list(), caseset, elig_controls, data = NULL, 
                       silent = FALSE) {
  entry <- eval(substitute(entry), data)
  exit <- eval(substitute(exit), data)
  fail <- eval(substitute(fail), data)
  origin <- eval(substitute(origin), data)
  elig_controls <- eval(substitute(elig_controls), data)
  n <- length(fail)
# individual id for identifying sets of split records
  ncc_id <- 1:n
# add this id to the supplied data frame for later merging
  data <- data.table(data)
  data[, "ncc_id" := ncc_id]
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
  
  grp <- rep(1, n)
  if (length(match) > 0) {
    grp <- as.numeric(factor(do.call(paste0, match)))
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
        n_eligible <- elig_controls[failure==TRUE]
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
# cases have probability of inclusion of 1 (comment this out, we still want
# the probability that they were ever included as a control)
#  tsplit[ncc_fail==TRUE, ncc_pr := 1]
  # if risk sets are smaller than requested number of cases, the calculated
# probability will be greater than 1. The correct value is 1.
  tsplit[ncc_pr > 1, ncc_pr := 1]
  # we no longer need the number of failures per set
  #tsplit[, nfail := NULL] 
  
  ncc_frame <- tsplit[, head(.SD, 1), by = list(ncc_id)]
  
  # coerce failure flag to numeric 
  # (logical var for case status might confuse people)
  ncc_frame[, ncc_fail := as.numeric(ncc_fail)]
  
  # merge ncc frame with original data
  setkey(ncc_frame, ncc_id)
  setkey(data, ncc_id)
  res <- data[ncc_frame]
#  res <- res[, c(names(ncc_frame), match_names), with=FALSE]
  res <- res[order(ncc_set,-ncc_fail), ]
  setkey(res, ncc_set)
  
  return(res)
}