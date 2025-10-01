check_taxonomy_consistency <- function(df) {
  library(dplyr)
  
  ranks <- c("k", "p", "c", "o", "f", "g", "s")
  inconsistencies <- list()
  
  for (i in seq_along(ranks)[-1]) {
    parent <- ranks[i - 1]
    child <- ranks[i]
    
    duplicates <- df %>%
      filter(!is.na(.data[[child]])) %>%
      distinct(.data[[child]], .data[[parent]]) %>%
      count(.data[[child]]) %>%
      filter(n > 1)
    
    if (nrow(duplicates) > 0) {
      inconsistencies[[child]] <- duplicates
    }
  }
  
  # Print results
  if (length(inconsistencies) > 0) {
    cat("⚠️ Inconsistencies found in the following taxonomic levels:\n\n")
    for (level in names(inconsistencies)) {
      cat(paste0("Level: ", level, "\n"))
      print(inconsistencies[[level]])
      cat("\n")
    }
  } else {
    cat("✅ No taxonomic inconsistencies found!\n")
  }
}
