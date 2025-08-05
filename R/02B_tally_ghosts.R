
tally_ghosts <- function(df){
hg_distinct <- df %>% select(sciname, ghost, n_seqs) %>% distinct()
num_ghosts <- sum(hg_distinct$ghost)
num_all <- nrow(hg_distinct)
pct_ghosts <- 100 * num_ghosts / num_all
mean_n_seqs <- mean(hg_distinct$n_seqs)
mean_n_seqs_non_ghosts <- mean(hg_distinct[hg_distinct$n_seqs > 0,]$n_seqs)
median_n_seqs_non_ghosts <- median(hg_distinct[hg_distinct$n_seqs > 0,]$n_seqs)

return(list(num_ghosts = num_ghosts,
            num_all = num_all,
            pct_ghosts = pct_ghosts,
            mean_n_seqs = mean_n_seqs,
            mean_n_seqs_non_ghosts = mean_n_seqs_non_ghosts,
            median_n_seqs_non_ghosts = median_n_seqs_non_ghosts))
}