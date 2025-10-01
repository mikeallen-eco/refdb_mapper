library(lme4)

fit_models_one_loo_type <- function(loo_method = "LOSO", 
                          assign_rubric = assign_rubric,
                          markers = markers,
                          loo_outcomes_list = outcomes){

  
all_pred_list <- lapply(1:length(markers), function (i){
  
# outcomes list for one marker
outcomes_marker <- loo_outcomes_list[markers[i]]
  
# create and empty list to hold all model predictions
pred_list <- list()

# loop through Incorrect, Abstain, and Correct metrics
if (loo_method %in% c("loso", "LOSO")) {
  metric_vector <- c("i", "a", "c")
  if (grepl(assign_rubric, pattern = "thresh")) {
    outcomes_df <- outcomes_marker[[1]]$loso
  }
  if (grepl(assign_rubric, pattern = "rdp")) {
    outcomes_df <- outcomes_marker[[1]]$loso_rdp
  }
  
} else{
  metric_vector <- c("i", "a")
  if (grepl(assign_rubric, pattern = "thresh")) {
    outcomes_df <- outcomes_marker[[1]]$lospo
  }
  if (grepl(assign_rubric, pattern = "rdp")) {
    outcomes_df <- outcomes_marker[[1]]$lospo_rdp
  }
  
}

for(j in metric_vector){  

# define model formula
  if (loo_method %in% c("loso", "LOSO")) {
    mod_formula <- as.formula(paste0(
      assign_rubric, "_", j,
      " ~ nnd + log(n_seqs) + nnd * log(n_seqs)"
    ))
  }
  
  if (loo_method %in% c("lospo", "LOSpO")) {
    mod_formula <- as.formula(paste0(
      assign_rubric, "_", j,
      " ~ nnd"
    ))
  }
  
  # fit model
mod_object <- glm(
  mod_formula,
  family = binomial,
  data = outcomes_df
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

pred_list[[j]] <- (list(pred_df = pred_df, mod = mod_object))

}

return(pred_list)
})

names(all_pred_list) <- markers

return(all_pred_list)

}
