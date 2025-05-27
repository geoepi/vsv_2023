get_tracer_stats <- function(log_file_path, use_burnin=0.10) {
  
  require(coda)
  
  mcmc_object <- get_log_file(log_file_path, burnin=use_burnin)
  
  if (!is.data.frame(mcmc_object)) {
    mcmc_object <- as.data.frame(as.mcmc(mcmc_object))
  }
  
  probs <- c(0.025, 0.5, 0.975)
  
  result_list <- list()
  
  numeric_columns <- sapply(mcmc_object, is.numeric)
  
  for(column_name in names(mcmc_object)[numeric_columns]) {
    quantiles <- quantile(mcmc_object[[column_name]], probs = probs, na.rm = TRUE)
    mean_value <- mean(mcmc_object[[column_name]], na.rm = TRUE)
    
    combined_stats <- c(quantiles, Mean = mean_value)
    
    result_list[[column_name]] <- combined_stats
  }
  
  result_df <- do.call(rbind, result_list)
  result_df <- apply(result_df, 2, function(x) round(x, 3))
  keep_rownames <- names(result_list)
  keep_colnames <- c("Q_0.025", "Q_0.5", "Q_0.975", "Mean")
  
  result_df <- as.data.frame(result_df)
  colnames(result_df) <- keep_colnames
  result_df$Parameter <- keep_rownames
  
  get_ess <- coda::effectiveSize(mcmc_object)
  result_df$ESS <- round(get_ess[rownames(result_df)], 0)
  
  result_df <- result_df %>%
    mutate(Median = Q_0.5) %>%
    select(Parameter, Mean, Median, Q_0.025, Q_0.975, ESS)
  
  rownames(result_df) <- NULL
  
  return(result_df)
}

# this part from https://github.com/laduplessis/beastio
get_log_file <- function(filename, burnin=use_burnin, maxsamples=-1, 
                         as.mcmc=TRUE, burninAsSamples=FALSE) {
  if (!burninAsSamples && burnin > 1) {
    stop("Error: Burnin must be a proportion of samples between 0 and 1.", call. = FALSE)
  }
  
  if (burninAsSamples && burnin != round(burnin)) {
    stop("Error: Burnin must be an integer number of states.", call. = FALSE)
  }
  
  logfile <- read.table(filename, sep="\t", header=TRUE, nrows=maxsamples)
  n <- nrow(logfile)
  
  if (!burninAsSamples) {
    burnSamples <- floor(burnin*n)
  } else {
    burnSamples <- burnin+1
  }
  
  if (burnSamples >= n) {
    stop("Error: Discarding all samples in the log file.", call. = FALSE)
  }
  
  logfile <- logfile[burnSamples:n,]
  
  if (as.mcmc == TRUE) {
    if (is.null(logfile$Sample)) {
      if (is.null(logfile$state)) {
        logfile$Sample <- as.numeric(rownames(logfile))
      } else {
        logfile$Sample <- logfile$state
        logfile$state  <- NULL
      }
    }
    
    for (i in names(logfile)) {
      if (all(is.na(logfile[[i]]))) logfile[[i]] <- NULL
    }
    
    start <- logfile$Sample[1]
    thin  <- logfile$Sample[2]-logfile$Sample[1]
    rownames(logfile) <- logfile$Sample
    logfile$Sample    <- NULL
    
    return(coda::mcmc(logfile, start=start, thin=thin))
    
  } else {
    return(logfile)
  }
}