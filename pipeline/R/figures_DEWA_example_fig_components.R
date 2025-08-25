source("mammals_Vences_16S.R")

# ---- example in one hydrobasin (Delaware River)

dewa <- -1529894232; rar <- -1529894602

# geotax <- fread("data/geotax.csv")
dewa_ghosts <- hydrobasin_ghosts %>%
  filter(HYBAS_ID %in% dewa)

(dewa_ghosts_tally <- tally_ghosts(dewa_ghosts))
dewa_info_summarized <- seq_info_summarized_by_hydrobasin %>%
  filter(HYBAS_ID %in% dewa)

dewa_hybas_pred_map_sf <- hybas_pred_map_sf %>%
  filter(HYBAS_ID %in% dewa)

#### MAP

library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)

# Filter your feature data
dewa_hybas_pred_map_sf <- hybas_pred_map_sf %>%
  filter(HYBAS_ID %in% dewa)

# Load world basemap from Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define orthographic projection centered on Eastern USA
ortho_proj <- "+proj=ortho +lat_0=40 +lon_0=-75 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Transform both the basemap and your spatial data
world_ortho <- st_transform(world, crs = ortho_proj)
dewa_ortho <- st_transform(dewa_hybas_pred_map_sf, crs = ortho_proj)
hybas_ortho <- st_transform(hybas_pred_map_sf, crs = ortho_proj)

# Create circular mask using a buffer around a center point
circle_mask <- st_point(c(0, 0)) |>
  st_sfc(crs = ortho_proj) |>
  st_buffer(dist = 6.4e6)  # ~6400 km = Earth's radius in meters

# Plot
ggplot() +
  geom_sf(data = circle_mask, fill = "gray90", color = NA) +  # background globe color
  geom_sf(data = world_ortho, fill = "gray50", color = "white", size = 0.2) +
  geom_sf(data = hybas_ortho, fill = "gray50", color = "gray80", size = 0.2) +
  geom_sf(data = dewa_ortho, fill = "red", alpha = 0.7, color = NA) +
  coord_sf(crs = ortho_proj) +
  theme_void() +
  theme(text = element_text(size = 18))

# Save figure
ggsave("figures/conceptual/globe_map_wo_hybas.png", height = 6, width = 6, dpi = 400)

### TABLE

dewa_hybas_error_data <- hybas_error_data %>%
  filter(HYBAS_ID %in% dewa)

fig_df <- dewa_hybas_error_data %>%
  select(order, geo_name, n_seqs, nnd, 
         p_misclass = preds_i, p_unclass = preds_a, p_correct = preds_c) %>%
  mutate(nnd = round(nnd),
         p_misclass = round(p_misclass, 1),
         p_unclass = round(p_unclass, 1),
         p_correct = round(p_correct, 1))

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
