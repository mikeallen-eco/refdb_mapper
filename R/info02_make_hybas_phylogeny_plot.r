library(ape)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
if (!exists("is.waive")) {
  is.waive <- function(x) inherits(x, "waiver")
}

make_hybas_phylogeny_plot <- function(hybas_data, # such as produced by get_polygon_attributes_from_coords()
                                 metric = "p_i",
                                 marker_list = markers[c(1,2,4,3,5)],
                                 tree_path = "data/phyl.tre"){
  
  info <- hybas_data
  
  phyl.tre <- read.tree(tree_path)
  
  # Keep only phyl_name + p_i columns
  traits <- info %>%
    mutate(
      mol_name = case_when(
        if_all(contains("n_seqs"), ~ .x == 0) ~ paste0("** ", mol_name, " **"),
        TRUE ~ mol_name
      )
    ) %>%
    select(phyl_name, mol_name, ends_with(metric)) %>%
    filter(!is.na(phyl_name))
  
  if(nrow(info)-nrow(traits) > 0){
  message("Note: removed ", nrow(info)-nrow(traits), " species that had no phylogeny name. These should be harmonized manually: ")
  message(paste(info[is.na(info$phyl_name),]$mol_name, collapse = "\n"))
  }
  
  # get duplicated phyl_name labels (MOL splits) to add into the phylogeny
  dup_tips <- traits %>%
    group_by(phyl_name) %>%
    filter(n() > 1)
  
  # Prune tree to just the tips in your info
  pruned_tre <- drop.tip(phyl.tre, setdiff(phyl.tre$tip.label, traits$phyl_name))
  
  # Start with pruned tree (keep all phyl_name tips)
  new_tree <- pruned_tre
  
  dup_tips <- traits %>%
    group_by(phyl_name) %>%
    filter(n() > 1)
  
# Loop over duplicates
for (phy_name in unique(dup_tips$phyl_name)) {
  
  # Get all mol_names for this phyl_name
  mols <- dup_tips %>% filter(phyl_name == !!phy_name) %>% pull(mol_name)
  
  # Skip the original tip name
  mols <- setdiff(mols, phy_name)
  
  # Find the tip number of the original phyl_name
  tip_id <- which(new_tree$tip.label == phy_name)
  
  # Add a zero-length tip for each extra mol_name
  for (mol in mols) {
    new_tree <- bind.tip(new_tree, tip.label = mol, edge.length = 0, where = tip_id)
  }
}

# Include new mol_name tips in traits_long
traits_expanded <- traits %>%
  bind_rows(
    dup_tips %>%
      mutate(phyl_name = mol_name) %>%   # assign new tip names
      select(names(traits))
  )

traits_long <- traits_expanded %>%
  pivot_longer(-c(phyl_name, mol_name), names_to = "marker", values_to = "metric")

# Define desired order of rings
desired_order <- paste0(marker_list, "_", metric)

# Convert marker column to factor with that order
traits_long <- traits_long %>%
  mutate(marker = factor(marker, levels = desired_order))

if(metric %in% c("p_i", "p_a")) {
  tip_colors <- traits_long %>%
    group_by(phyl_name) %>%
    summarise(tip_metric = min(metric, na.rm = TRUE))
} else{
  tip_colors <- traits_long %>%
    group_by(phyl_name) %>%
    summarise(tip_metric = max(metric, na.rm = TRUE))
}

traits_long <- traits_long %>%
  mutate(metric = case_when(metric %in% 0 ~ NA,
                            TRUE ~ metric))

# Create a lookup table for tip labels
tip_label_lookup <- traits_expanded %>%
  select(phyl_name, mol_name) %>%
  distinct() %>%
  # For tips where mol_name == phyl_name, just keep the name
  mutate(label_to_use = mol_name)

# Join with tip_colors so we can map both color and label
tip_colors_labels <- tip_colors %>%
  left_join(tip_label_lookup, by = "phyl_name")

# Compute global min/max across all metric values (for tips and heatmap)
p_min <- min(traits_long$metric, tip_colors$tip_metric, na.rm = TRUE)
p_max <- max(traits_long$metric, tip_colors$tip_metric, na.rm = TRUE)

# Start base plot
p <- suppressWarnings(
  ggtree(new_tree, layout = "circular") %<+% tip_colors_labels
)

# Make terminal branches longer
p <- p + hexpand(10)

# Number of markers / rings
n_rings <- length(unique(traits_long$marker))

# Estimate total radial extension of heatmap rings
# geom_fruit width units are roughly proportional to tree height; add a generous margin
total_ring_offset <- 0.05 + n_rings * 6 * 0.5  + 19 # previous offset + half of total widths

# Add tip labels using mol_name
p <- p + geom_tiplab(
  aes(label = label_to_use, color = tip_metric),
  size = 3,
  align = TRUE,
  offset = total_ring_offset + 10
)

# Add heatmap rings
p <- p + geom_fruit(
  data = traits_long,
  geom = geom_tile,
  mapping = aes(y = phyl_name, x = marker, fill = metric),
  width = 6,
  offset = 0
)

# create legend title
if(metric %in% "p_i"){legend_title <- "p(misclass.)"}
if(metric %in% "p_c"){legend_title <- "p(correct)"}
if(metric %in% "p_a"){legend_title <- "p(unclass.)"}

# Use a single unified scale for both tip color and heatmap fill
shared_scale <- scale_fill_viridis_c(
  option = "C",
  name = legend_title,
  limits = c(p_min, p_max)
)

p <- p +
  shared_scale +
  scale_color_viridis_c(
    option = "C",
    limits = c(p_min, p_max),
    guide = "none"  # hide duplicate legend
  )
# p

# Stack the ring labels in the middle, spaced vertically
ring_labels_center <- data.frame(
  label = rev(marker_list),        # inside → outside order
  y = seq_along(marker_list)
)

p <- p +
  annotate(
    "text",
    x = 0, y = 0,                   # center of the plot
    label = paste(rev(marker_list), collapse = "\n"),
    size = 4, hjust = 0.5, vjust = 0.5
  )

# Stack the ring labels in the middle, spaced vertically
# ring_labels_center <- data.frame(
#   label = rev(marker_list),        # inside → outside order
#   y = seq_along(marker_list)
# )
# 
# p <- p + 
#   geom_label(
#     data = ring_labels_center,
#     aes(x = 0, y = y, label = label),
#     inherit.aes = FALSE,
#     size = 4,
#     hjust = 0.5, vjust = 0.5,
#     fill = "white", alpha = 0.7,   # semi-transparent white background
#     label.size = NA                # removes the border box
#   )


return(p)

}
