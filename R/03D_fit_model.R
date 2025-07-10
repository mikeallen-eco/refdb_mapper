library(lme4)

fit_model <- function(assign_method, loo_method, acc_metric, loo_outcomes_df){

# define model formula
  if (loo_method %in% "LOSO") {
    mod_formula <- as.formula(paste0(
      acc_metric,
      " ~ nnd + log(n_seqs) + nnd * log(n_seqs)"
    ))
  }
  
  if (loo_method %in% "LOSpO") {
    mod_formula <- as.formula(paste0(
      acc_metric,
      " ~ nnd"
    ))
  }
  
  # fit model
mod_object <- glm(
  mod_formula,
  family = binomial,
  data = loo_outcomes_df %>% filter(method %in% assign_method)
); summary(mod_object)

preds <- predict(
  mod_object,
  newdata = expand.grid(
    nnd = seq(0, 30, length.out = 100),
    n_seqs = seq(2, 11, length.out = 100)
  ),
  type = "response",
  re.form = NA
)

pred_df <- expand.grid(nnd = seq(0, 30, length.out = 100),
                       n_seqs = seq(2, 11, length.out = 100)) %>%
  mutate(preds = preds)

return(list(pred_df = pred_df,
            mod = mod_object))

}
