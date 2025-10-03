mia_metrics <- function(phyobj) {
  # ensure required packages are loaded
  if (!requireNamespace("mia", quietly = TRUE)) {
    stop("Package 'mia' is required but not installed.")
  }
  
  rc <- readcount(phyobj)
  ab <- abundances(phyobj)
  n  <- ntaxa(phyobj)
  
  safe_ntaxa <- function(expr) {
    tryCatch(ntaxa(expr), error = function(e) NA)
  }
  
  # collect metrics in a named list
  metrics <- list(
    min_read      = min(rc),
    max_read      = max(rc),
    total_reads   = sum(rc),
    mean_reads    = round(mean(rc), 0),
    median_reads  = median(rc),
    total_asvs    = n,
    singletons    = safe_ntaxa(rare(phyobj, detection = 1, prevalence = 0)),
    singletons_pc = tryCatch(
      round(100 * safe_ntaxa(rare(phyobj, detection = 1, prevalence = 0)) / n, 3),
      error = function(e) NA
    ),
    sparsity      = round(mean(ab == 0), 3)
  )
  
  return(as.data.frame(metrics))
}