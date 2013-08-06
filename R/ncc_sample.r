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
#' @param silent if \code{FALSE}, echoes a fullstop to the screen as each risk 
#' set is generated, and as each case-control set is sampled. If TRUE,
#' suppresses output to the terminal.
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
#' Control selection is performed using the \code{\link{sample}} function. It is
#' therefore important to set the seed for the pseudo random number generator 
#' (\code{\link{set.seed}}) to ensure reproducibility.
#' 
#' @return
#' A \code{data.frame} comprising:
#' \item{ncc_id}{a unique individual identifier}
#' \item{ncc_set}{a case-control set identifier}
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
#' @export
#' 
ncc_sample <- function(entry = 0, exit, fail, origin = 0, controls = 1, 
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
data$ncc_id <- ncc_id
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
  set <- 0
  nomatch <- 0
  incomplete <- 0
  ties <- 0 
  fg <- unique(grp[fail != 0])
  gresult <- vector("list", length(fg))
  # index for groups
  gii <- 0
  # iterator for displaying dots
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
        df <- cbind(ncc_id, 
                    set, 
                    failure, 
                    rep(sum(failure), length(failure)),
                    rep(n_eligible, length(failure)),
                    rep(t_end, length(failure)))
        colnames(df) <- c("ncc_id", "ncc_set", "ncc_fail", 
                          "nfail", "ncc_elig_co", "ncc_time")
#         df <- data.frame(ncc_id   = ncc_id, 
#                          ncc_set  = set, 
#                          ncc_fail = failure, 
#                          nfail    = rep(sum(failure), length(failure)),
#                          ncc_elig_co = rep(n_eligible, length(failure)), 
#                          ncc_time     = rep(t_end, length(failure)))
        sresult[[tjj]] <- df[!ineligible, ]
      }
      tresult[[tii]] <- do.call("rbind", sresult)
    }
    gresult[[gii]] <- do.call("rbind", tresult)
    rm(tresult)
  }
  tsplit <- as.data.frame(do.call("rbind", gresult))
# Calculate probability of inclusion 
# (Samuelsen, Biometrika, 1997, 84(2) p379)
  pr_not <- 1 - ((tsplit$nfail * controls) / tsplit$ncc_elig_co)
  product <- aggregate(pr_not, by = list(tsplit$ncc_id), prod)
  names(product) <- c("ncc_id", "prod_pr_not")
  pr <- data.frame(product["ncc_id"], product["prod_pr_not"])
  pr$ncc_pr <- 1 - pr$prod_pr_not
  pr$prod_pr_not <- NULL
  tsplit <- merge(tsplit, pr, sort = FALSE)
# cases have probability of inclusion of 1
  tsplit$ncc_pr[tsplit$ncc_fail==TRUE] <- 1
# if risk sets are smaller than requested number of cases, the calculated
# probability will be greater than 1. The correct value is 1.
  tsplit$ncc_pr[tsplit$ncc_pr > 1] <- 1
  # we no longer need the number of failures per set
  tsplit$nfail <- NULL
  
# sample controls from risk sets  
  if (!silent) {
    cat("\nSampling from risk sets:\n")
  }
  nsets <- max(tsplit$ncc_set)
  samp <- vector("list", nsets)
  for (s in 1:nsets) {
  if (!silent) {
    dots <- ifelse((s %% 50 == 0), ". 50\n", ".")
    cat(dots)
  }
    riskset <-tsplit[tsplit$ncc_set==s, ]
    ncase <- sum(riskset$ncc_fail)
    nelig <- sum(!riskset$ncc_fail)
    ncont <- controls * ncase
    if (nelig < ncont) {
      ncont <- nelig
      if (ncont > 0) {
        incomplete <- incomplete + 1
      }
    }
    if (ncont == 0) {
      nomatch <- nomatch + ncase
    }
    else {
      co_samp <- sample((1:nrow(riskset))[!riskset$ncc_fail], ncont)
      ca_samp <- (1:nrow(riskset))[riskset$ncc_fail]
      samp_index <- c(ca_samp, co_samp)
      samp[[s]] <- riskset[samp_index, ]
    }
  }
  ncc_frame <- do.call("rbind.data.frame", samp)  
  # coerce failure flag to numeric 
  # (logical var for case status might confuse people)
  ncc_frame$ncc_fail <- as.numeric(ncc_frame$ncc_fail)
  res <- merge(ncc_frame, data, by="ncc_id")
  res <- res[, c(names(ncc_frame), match_names, inc_names)]
  res <- res[order(res$ncc_set,res$ncc_fail), ]
  row.names(res) <- 1:nrow(res)
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
  return(data.frame(res))
}
  