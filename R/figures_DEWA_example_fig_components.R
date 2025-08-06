source("mammals_Vences_16S.R")

# ---- Step 5 example in one hydrobasin (Delaware River)

dewa <- -1529894232; rar <- -1529894602

# geotax <- fread("data/geotax.csv")
dewa_ghosts <- hydrobasin_ghosts %>%
  filter(HYBAS_ID %in% dewa)

(dewa_ghosts_tally <- tally_ghosts(dewa_ghosts))
dewa_info_summarized <- seq_info_summarized_by_hydrobasin %>%
  filter(HYBAS_ID %in% dewa)

dewa_hybas_pred_map_sf <- hybas_pred_map_sf %>%
  filter(HYBAS_ID %in% dewa)

# hybas_NE <- subset_hydrobasins_states(hydrobasin_map = hybas_pred_map_sf)

# Get US states from Natural Earth
us_states <- ne_states(country = "United States of America", returnclass = "sf")

# List of Northeastern states (per US Census Bureau definition)
northeast_states <- c(
  "Connecticut", "Maine", "Massachusetts", "New Hampshire", "Virginia", "West Virginia",
  "Rhode Island", "Vermont", "New Jersey", "New York", "Pennsylvania", "Maryland", "Delaware"
)

# Filter for northeastern states
ne_us_sf <- us_states %>%
  filter(name %in% northeast_states)

ggplot(ne_us_sf) +
  geom_sf() +
  geom_sf(data = dewa_hybas_pred_map_sf, fill = "red", alpha = 0.7) +
  theme_minimal() +
  theme(text = element_text(size = 18))

ggsave("figures/conceptual/NE_map.png", height = 6, width = 6, dpi = 400)

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
