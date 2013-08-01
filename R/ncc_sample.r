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
  nn <- (1 + controls) * sum(fail != 0)
  pr <- numeric(nn)
  sr <- numeric(nn)
  tr <- vector("numeric", nn)
  fr <- numeric(nn)
  nn <- 0
  if (!silent) {
    cat("\nGenerating risk sets: ")
  }
  set <- 0
  nomatch <- 0
  incomplete <- 0
  ties <- 0 
  fg <- unique(grp[fail != 0])
  gresult <- vector("list", length(fg))
  # index for groups
  gii <- 0
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
        if (!silent) {
          cat(".")
        }
        set <- set + 1
        failure <- case & case_index==tjj
        noncase <- (grp == g) & (t_entry <= t_end) & (t_exit >= t_end) & !failure
        n_eligible <- sum(noncase)
        ineligible <- !failure & !noncase
        df <- data.frame(ncc_id   = ncc_id, 
                         ncc_set  = set, 
                         ncc_fail = failure, 
                         nfail    = rep(sum(failure), length(failure)),
                         ncc_elig_co = rep(n_eligible, length(failure)), 
                         ncc_time     = rep(t_end, length(failure)))
        sresult[[tjj]] <- df[!ineligible, ]
      }
      tresult[[tii]] <- do.call(rbind.data.frame, sresult)
    }
    gresult[[gii]] <- do.call(rbind.data.frame, tresult)
    rm(tresult)
  }
  tsplit <- do.call(rbind.data.frame, gresult)
# Calculate probability of inclusion 
# (Samuelsen, Biometrika, 1997, 84(2) p379)
  pr_not <- 1 - ((tsplit$nfail * controls) / tsplit$ncc_elig_co)
  product <- aggregate(pr_not, by = list(tsplit$ncc_id), prod)
  names(product) <- c("ncc_id", "prod_pr_not")
  pr <- data.frame(product["ncc_id"], product["prod_pr_not"])
  pr$ncc_pr <- 1 - pr$prod_pr_not
  pr$prod_pr_not <- NULL
  tsplit <- merge(tsplit, pr, sort = FALSE)
  tsplit$ncc_pr[tsplit$ncc_fail==TRUE] <- 1
  # we no longer need the number of failures per set
  tsplit$nfail <- NULL
  
  
# sample controls from risk sets  
  nsets <- max(tsplit$ncc_set)
  samp <- vector("list", nsets)
  for (s in 1:nsets) {
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
  ncc_frame <- do.call(rbind.data.frame, samp)  
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
  