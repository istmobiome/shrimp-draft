summarize_objs <- function(objs) {
  results <- data.frame(
    me_dataset    = character(),
    total_asvs    = integer(),
    total_reads   = numeric(),
    total_samples = integer(),
    stringsAsFactors = FALSE
  )
  
  for (obj in objs) {
    if (exists(obj)) {
      x <- get(obj)
      
      # safely extract values
      asvs <- tryCatch(nrow(x$tax_table), error = function(e) NA)
      reads <- tryCatch(sum(x$sample_sums()), error = function(e) NA)
      samples <- tryCatch(nrow(x$sample_table), error = function(e) NA)
      
      if (!is.na(asvs)) {
        results <- rbind(
          results,
          data.frame(
            me_dataset = obj,
            total_asvs = asvs,
            total_reads = reads,
            total_samples = samples,
            stringsAsFactors = FALSE
          )
        )
      } else {
        message(sprintf("Skipping '%s' (missing or invalid tax_table)", obj))
      }
    } else {
      message(sprintf("Skipping '%s' (object missing)", obj))
    }
  }
  
  results
}