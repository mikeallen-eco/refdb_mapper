# data vis & info by hydrobasin coordinates
source("R/setup_present.R")

# ---- Step 1: get marker detectability data for a hydrobasin from coordinates

# nj_data <- get_polygon_attributes_from_coords(-74.759831, 40.568704)
ct_data <- get_polygon_attributes_from_coords(-72.92116, 41.31607, include_marine = F)
pf_data <- get_polygon_attributes_from_coords(-59.912786, -2.069209, include_marine = F)

# ---- Step 2: choose best combinations of markers

ct_best <- pick_best_markers(ct_data)
pf_best <- pick_best_markers(pf_data)

plot_best_markers(ct_best, metric = "c")
ggsave("figures/best_markers_p_c_CT.png", height = 10, width = 10, dpi = 400)
plot_best_markers(pf_best)
ggsave("figures/best_markers_p_c_PF.png", height = 10, width = 10, dpi = 400)

plot_best_markers(ct_best, metric = "i")
ggsave("figures/best_markers_p_i_CT.png", height = 10, width = 10, dpi = 400)
plot_best_markers(pf_best, metric = "i")
ggsave("figures/best_markers_p_i_PF.png", height = 10, width = 10, dpi = 400)

# ---- Step 3: plot % correct vs. % incorrect

scatterplot_best_markers(ct_best, top_n = 10, num_markers = 1:2)
ggsave("figures/best_scatter_CT_top_n_all.png", height = 6, width = 7, dpi = 400)

scatterplot_best_markers(pf_best, top_n = 10, num_markers = 1:2)
ggsave("figures/best_scatter_PF_top_n_10.png", height = 6, width = 7, dpi = 400)

scatterplot_best_markers(ct_best, top_n = 5, num_markers = 1:2)
ggsave("figures/best_scatter_CT_top_n_5.png", height = 6, width = 7, dpi = 400)

scatterplot_best_markers(pf_best, top_n = 5, num_markers = 1:2)
ggsave("figures/best_scatter_PF_top_n_5.png", height = 6, width = 7, dpi = 400)

# ---- Step 4: phylogenetic visualization of detectability by marker

make_hybas_phylogeny_plot(hybas_data = ct_data, metric = "p_c", rubric = "rdp90")
ggsave("figures/circular_phylogeny_rdp90_p_c_CT.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = ct_data, metric = "p_c", rubric = "blast98")
ggsave("figures/circular_phylogeny_blast98_p_c_CT.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = pf_data, metric = "p_c", rubric = "rdp90")
ggsave("figures/circular_phylogeny_rdp90_p_c_PF.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = pf_data, metric = "p_c", rubric = "blast98")
ggsave("figures/circular_phylogeny_blast98_p_c_PF.png", height = 12, width = 12, dpi = 400)

make_hybas_phylogeny_plot(hybas_data = ct_data, metric = "p_i", rubric = "rdp90")
ggsave("figures/circular_phylogeny_rdp90_p_i_CT.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = ct_data, metric = "p_i", rubric = "blast98")
ggsave("figures/circular_phylogeny_blast98_p_i_CT.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = pf_data, metric = "p_i", rubric = "rdp90")
ggsave("figures/circular_phylogeny_rdp90_p_i_PF.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = pf_data, metric = "p_i", rubric = "blast98")
ggsave("figures/circular_phylogeny_blast98_p_i_PF.png", height = 12, width = 12, dpi = 400)

make_hybas_phylogeny_plot(hybas_data = ct_data, metric = "p_a", rubric = "rdp90")
ggsave("figures/circular_phylogeny_rdp90_p_a_CT.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = nj_data, metric = "p_a", rubric = "blast98")
ggsave("figures/circular_phylogeny_blast98_p_a_CT.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = pf_data, metric = "p_a", rubric = "rdp90")
ggsave("figures/circular_phylogeny_rdp90_p_a_PF.png", height = 12, width = 12, dpi = 400)
make_hybas_phylogeny_plot(hybas_data = pf_data, metric = "p_a", rubric = "blast98")
ggsave("figures/circular_phylogeny_blast98_p_a_PF.png", height = 12, width = 12, dpi = 400)





# ---- Step 0: plot error rate model effects by marker

fits <- readRDS("data/fits_20251003.rds")
eplots <- plot_predicted_loso_lopso_error(preds = fits, markers = markers, rubrics = rubrics)
eplots$blast98$RiazVert1_12S$loso$i
eplots$blast98$Vences_16S$loso$i
eplots$rdp90$RiazVert1_12S$loso$i
eplots$rdp90$Vences_16S$loso$i

eplots$blast98$RiazVert1_12S$loso$a
eplots$blast98$Vences_16S$loso$a
eplots$rdp90$RiazVert1_12S$loso$a
eplots$rdp90$Vences_16S$loso$a

eplots$blast98$RiazVert1_12S$loso$c
eplots$blast98$Vences_16S$loso$c
eplots$rdp90$RiazVert1_12S$loso$c
eplots$rdp90$Vences_16S$loso$c


eplots$blast98$RiazVert1_12S$lospo$i
eplots$blast98$Vences_16S$lospo$i
eplots$rdp90$RiazVert1_12S$lospo$i
eplots$rdp90$Vences_16S$lospo$i


eplots$blast98$RiazVert1_12S$lospo$a
eplots$blast98$Vences_16S$lospo$a
eplots$rdp90$RiazVert1_12S$lospo$a
eplots$rdp90$Vences_16S$lospo$a

### TABLE

hid_hybas_error_data <- nj_data

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


# ---- Step x: global locator map for hydrobasin coordinates

make_global_locator_map(-72.92116, 41.31607) # Peabody Museum, CT
ggsave("figures/globe_locator_map_CT.png", height = 6, width = 6, dpi = 400)
make_global_locator_map(lon = -59.912786, lat = -2.069209)
ggsave("figures/globe_locator_map_PF.png", height = 6, width = 6, dpi = 400)
