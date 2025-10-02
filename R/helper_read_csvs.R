# function to compile all the csvs in a folder into one dataframe

library(dplyr)

read_csvs <- function(directory,
                      pattern_str = "*.csv") {
  filenames <- list.files(directory, pattern = pattern_str, full.names = T)
  combined_list <- lapply(seq_along(filenames),
                          function(x){
                           df <- read.csv(filenames[x])
                           if("acc" %in% colnames(df)){df$acc <- as.character(df$acc)} # patch for rare cases where accession number is just a number
                           return(df)
                          })

  combined_df <- do.call(bind_rows, combined_list)
  return(combined_df)
}
