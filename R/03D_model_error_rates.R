model_error_rates <- function(loso_GB_compiled){

loso_outcomes <- compile_loso_ghostblaster_outcomes(loso_GB_compiled)$b %>%
  mutate(method = "BLAST") %>%
  bind_rows(compile_loso_ghostblaster_outcomes(loso_GB_compiled)$gb %>% 
              mutate(method = "GhostBLASTer")) %>%
  select(qseqid, method, starts_with("top_match"), starts_with("score"), starts_with("max"), true_ncbi_name:conf90_a)

# prepare BLAST data for binomial model
df <- loso_outcomes %>%
  dplyr::rename(n_seqs = n) %>%
  group_by(method, true_ncbi_name, order, nnd, n_seqs) %>%
  summarize(thresh98 = sum(thresh98_i),
            ecotag_i = sum(ecotag_i),
            conf90_i = sum(conf90_i),
            thresh98_i = sum(thresh98_i),
            ecotag_c = sum(ecotag_c),
            thresh99_i = sum(thresh99_i),
            thresh98_c = sum(thresh98_c),
            conf90_c = sum(conf90_c),
            n = length(true_ncbi_name),
            .groups = "drop") %>%
  mutate(nndS = wiqid::standardize(nnd),
         n_seqsS = wiqid::standardize(n_seqs))

library(lme4)

m_thresh99_i <- glm(
  cbind(thresh99_i, n - thresh98_i) ~ nnd + log(n_seqs) + nnd*log(n_seqs),
  family = binomial,
  data = df %>% filter(method %in% "BLAST")
); summary(m_thresh99_i)

m_ecotag_i <- glm(
  cbind(ecotag_i, n - thresh98_i) ~ nnd + log(n_seqs) + nnd*log(n_seqs),
  family = binomial,
  data = df %>% filter(method %in% "BLAST")
); summary(m_ecotag_i)

m_conf90_i <- glm(
  cbind(conf90_i, n - conf90_i) ~ nnd + log(n_seqs) + nnd*log(n_seqs),
  family = binomial,
  data = df %>% filter(method %in% "BLAST")
); summary(m_conf90_i)

m_conf90_gb_i <- glm(
  cbind(conf90_i, n - conf90_i) ~ nnd + log(n_seqs) + nnd*log(n_seqs),
  family = binomial,
  data = df %>% filter(method %in% "GhostBLASTer")
); summary(m_conf90_gb_i)

preds_thresh99_i <- predict(
  m_thresh99_i,
  newdata = expand.grid(
    nnd = seq(0, 30, length.out = 100),
    n_seqs = seq(2, 11, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_ecotag_i <- predict(
  m_ecotag_i,
  newdata = expand.grid(
    nnd = seq(0, 30, length.out = 100),
    n_seqs = seq(2, 11, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_conf90_i <- predict(
  m_conf90_i,
  newdata = expand.grid(
    nnd = seq(0, 30, length.out = 100),
    n_seqs = seq(2, 11, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

preds_conf90_gb_i <- predict(
  m_conf90_gb_i,
  newdata = expand.grid(
    nnd = seq(0, 30, length.out = 100),
    n_seqs = seq(2, 11, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

pred_df_thresh99_i <- expand.grid(nnd = seq(0, 30, length.out = 100),
                                  n_seqs = seq(2, 11, length.out = 100)) %>%
  mutate(preds = preds_thresh99_i)

pred_df_ecotag_i <- expand.grid(nnd = seq(0, 30, length.out = 100),
                                n_seqs = seq(2, 11, length.out = 100)) %>%
  mutate(preds = preds_ecotag_i)

pred_df_conf90_i <- expand.grid(nnd = seq(0, 30, length.out = 100),
                                n_seqs = seq(2, 11, length.out = 100)) %>%
  mutate(preds = preds_conf90_i)

pred_df_conf90_gb_i <- expand.grid(nnd = seq(0, 30, length.out = 100),
                                n_seqs = seq(2, 11, length.out = 100)) %>%
  mutate(preds = preds_conf90_gb_i)


return(list(pred_df_thresh99_i = pred_df_thresh99_i,
            pred_df_ecotag_i = pred_df_ecotag_i,
            pred_df_conf90_i = pred_df_conf90_i,
            pred_df_conf90_gb_i = pred_df_conf90_gb_i,
            mod_thresh99_i = m_thresh99_i,
            mod_ecotag_i = m_ecotag_i,
            mod_conf90_i = m_conf90_i,
            mod_conf90_gb_i = m_conf90_gb_i))
}
