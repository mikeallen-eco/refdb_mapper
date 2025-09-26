# data vis & info by hydrobasin coordinates
source("R/setup_info.R")

# ---- Step 1: global locator map for hydrobasin coordinates

make_global_locator_map(-74.759831, 40.568704) # Readington, NJ -59.912786, -2.069209 = Presidente Figueiredo, BR
ggsave("figures/globe_locator_map_raritan.png", height = 6, width = 6, dpi = 400)
make_global_locator_map(lon = -59.912786, lat = -2.069209)
ggsave("figures/globe_locator_map_raritan_PF.png", height = 6, width = 6, dpi = 400)


# ---- Step 3: phylogenetic vis of attributes
# note: make them work on "sf_data" created using get_attributes... function once

make_hybas_phylogeny_plot(metric = "p_c")
ggsave("figures/circular_phylogeny_test_p_c_3.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(lon = -59.912786, lat = -2.069209, metric = "p_c")
ggsave("figures/circular_phylogeny_test_p_c_PF.png", height = 12, width = 12, dpi = 400)

make_hybas_phylogeny_plot(metric = "p_i")
ggsave("figures/circular_phylogeny_test_p_i_1.png", height = 12, width = 12, dpi = 400)

make_hybas_phylogeny_plot(metric = "p_a")
ggsave("figures/circular_phylogeny_test_p_a_2.png", height = 12, width = 12, dpi = 400)

info <- get_polygon_attributes_from_coords(-74.759831, 40.568704, polygons = final_sf)

head(info)

# Install if needed
# BiocManager::install("ggtreeExtra")

library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(dplyr)
library(tidyr)

phyl.tre <- read.tree("data/phyl.tre")

# Keep only phyl_name + p_i columns
traits <- info %>%
  select(mol_name, phyl_name, ends_with("p_i")) %>%
  filter(!is.na(phyl_name))

# Prune tree to just the tips in your info
pruned_tre <- drop.tip(phyl.tre, setdiff(phyl.tre$tip.label, traits$phyl_name))

# Reshape to long format for heatmap plotting
traits_long <- traits %>%
  pivot_longer(-phyl_name, names_to = "marker", values_to = "p_i")

# Start base plot
p <- suppressWarnings(ggtree(pruned_tre, layout = "circular"))

# Add heatmap rings (all 5 in one layer)
p <- p + geom_fruit(
  data = traits_long,
  geom = geom_tile,
  mapping = aes(y = phyl_name, x = marker, fill = p_i),
  width = 5,        # thickness of each ring
  offset = 0.1        # distance away from tree
)

# Continuous scale
p <- p + scale_fill_viridis_c(option = "C", name = "p_i values")

p



### TABLE

hid_hybas_error_data <- final_sf %>%
  filter(HYBAS_ID %in% hid)

fig_df_n_seqs <- hid_hybas_error_data %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  select(order, mol_name, contains("n_seqs"))

fig_df_p_i <- hid_hybas_error_data %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  select(order, mol_name, nnd,
         contains("p_i")) %>%
  mutate(nnd = round(nnd)) %>%
  mutate(across(contains("p_i"), ~ round(.*100, 1)))


library(knitr)
library(kableExtra)

# Optional: Rename columns to be more readable
fig_df_pretty <- fig_df %>%
  arrange(desc(p_misclass)) %>%
  dplyr::rename(
    `Order` = order,
    `Species` = geo_name,
    `# Sequences` = n_seqs,
    `NND` = nnd,
    `P(Misclass)` = p_misclass,
    `P(Unclass)` = p_unclass,
    `P(Correct)` = p_correct
  )

fig_df_pretty_least <- fig_df_pretty %>%
  arrange(`P(Misclass)`)

fig_df_pretty_alpha <- fig_df_pretty %>%
  arrange(`Species`)

# Print nice table
# most error prone
most_table_html <- knitr::kable(fig_df_pretty, format = "html", digits = 3, escape = TRUE) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

# least error prone
least_table_html <- knitr::kable(fig_df_pretty_least, format = "html", digits = 3, escape = TRUE) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

# alphabetical
alpha_table_html <- knitr::kable(fig_df_pretty_alpha, format = "html", digits = 3, escape = TRUE) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = FALSE)

# Save to file
kableExtra::save_kable(most_table_html, "figures/conceptual/eg_most_error_table.html")
kableExtra::save_kable(least_table_html, "figures/conceptual/eg_least_error_table.html")
kableExtra::save_kable(alpha_table_html, "figures/conceptual/eg_alpha_table.html")
