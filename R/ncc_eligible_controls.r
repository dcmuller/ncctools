#' @title Calculate the number of eligible controls for each case for a nested 
#' case control study
#' 
#' @description
#' \code{ncc_eligible_controls} calculates the number of eligible controls for each
#' case for a nested case control study. Given time of entry, time of exit, and 
#' exit status, 
#' risk sets are computed at each failure time. If matching variables are specified, 
#' \code{ncc_eligible_controls} creates matched or stratified risk sets, 
#' which are computed separately within matching strata. 
#' 
#' @param entry time of entry to follow-up
#' @param exit time of exit from follow-up
#' @param fail indicator of status on exit from follow-up (censored=0, fail=1)
#' @param origin the origin of the analysis time-scale. For instance, date of
#' birth, if age is the desired time-scale.
#' @param match a list of categorical variables for matching cases and controls
#' @param data a \code{\link{data.frame}} which contains the follow-up and
#' matching variables.
#' @param include a list of variables to include in the returned data.frame
#' @param silent if \code{FALSE}, provides entertainment by echoing a fullstop 
#' to the screen as each risk set is generated. If TRUE, output to the
#' console is suppressed.
#' 
#' @details
#' Given follow-up information from a cohort study, \code{ncc_eligible controls}
#'  generates
#' risk sets at each observed failure time, and counts the number of eligible 
#' controls within each risk set. This information can be used to calculate the
#' probability of inclusion in the NCC study (see \code{\link{ncc_selection_prob}})
#' 
#' 
#' @return
#' A \code{data.frame} comprising:
#' \item{ncc_id}{a unique individual identifier}
#' \item{ncc_elig_co}{a count of the number of controls eligible for selection
#' in the set}
#' 
#' followed by the variables specified in the \code{match} and \code{include} 
#' lists.
#' 
#' @author David C Muller
#' 
#' @import data.table
#' @export
#' 
ncc_eligible_controls <- function(entry = 0, exit, fail, origin = 0, 
                       match = list(), include = list(), data = NULL, 
                       silent = FALSE) {
  entry <- eval(substitute(entry), data)
  exit <- eval(substitute(exit), data)
  fail <- eval(substitute(fail), data)
  origin <- eval(substitute(origin), data)
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
        sresult[[tjj]]<- df[df$ncc_fail==TRUE,,drop=FALSE]
      }
      tresult[[tii]] <- rbindlist(sresult)
    }
    gresult[[gii]] <- rbindlist(tresult)
    rm(tresult)
  }
  tsplit <- rbindlist(gresult)
  
  # we no longer need the number of failures per set
  #tsplit[, nfail := NULL] 
  
  ncc_frame <- tsplit
  
  # coerce failure flag to numeric 
  # (logical var for case status might confuse people)
  ncc_frame[, ncc_fail := as.numeric(ncc_fail)]

  
  # merge ncc frame with original data
  setkey(ncc_frame, ncc_id)
  setkey(data, ncc_id)
  res <- data[ncc_frame]
  res <- res[, c(names(ncc_frame), match_names, inc_names), with=FALSE]
  res <- res[order(ncc_set,-ncc_fail), ]
  setkey(res, ncc_set)
  
  return(res)
}
  