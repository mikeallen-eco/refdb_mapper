library(purrr)

# reorganize error model list structure to match: reorg$marker1$loso$i, reorg$marker1$lospo$i
reorganize_preds <- function(preds_loso, preds_lospo) {
  # get union of marker names across both
  markers <- union(names(preds_loso), names(preds_lospo))
  
  # build reorganized list
  map(markers, function(m) {
    list(
      loso  = preds_loso[[m]],
      lospo = preds_lospo[[m]]
    )
  }) |> set_names(markers)
}