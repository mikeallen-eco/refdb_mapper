library(shiny)
library(sf)
library(leaflet)
library(dplyr)
library(viridis)
library(DT)

# -------------------------------
# Load data
# -------------------------------
species_data <- readRDS("data/hydrobasin_refdb_nnd_info.rds")
hydrobasins <- readRDS("data/hydrobasin_map_simple.rds")

# Keep only polygons/multipolygons
hydrobasins <- hydrobasins %>%
  filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON"))

# Ensure proper geometry
st_geometry(hydrobasins) <- st_geometry(hydrobasins)

# -------------------------------
# Precompute stats per hydrobasin
# -------------------------------
ghost_stats <- species_data %>%
  group_by(HYBAS_ID) %>%
  summarize(
    n_species = n_distinct(sciname),
    n_ghosts = sum(V5_12S == 0 & MiMamm_12S == 0 & Vences_16S == 0),
    pct_ghosts = 100 * n_ghosts / n_species,
    .groups = "drop"
  )

hydrobasins <- hydrobasins %>%
  left_join(ghost_stats, by = "HYBAS_ID")

# -------------------------------
# Color palette
# -------------------------------
pal <- colorNumeric(
  palette = "inferno",
  domain = hydrobasins$pct_ghosts,
  na.color = "transparent"
)

# -------------------------------
# UI
# -------------------------------
ui <- fluidPage(
  titlePanel("Reference Database Explorer"),
  fluidRow(
    column(7, leafletOutput("map", height = "80vh")),
    column(5, DTOutput("basin_table"))
  )
)

# -------------------------------
# Server
# -------------------------------
server <- function(input, output, session) {
  
  # Render the leaflet map
  output$map <- renderLeaflet({
    leaflet(hydrobasins) %>%
      addProviderTiles("CartoDB.Positron") %>%
      addPolygons(
        layerId = ~HYBAS_ID,
        weight = 0.5,          # thinner borders
        color = "#555",        # slightly darker but transparent enough
        fillColor = ~pal(pct_ghosts),
        fillOpacity = 0.6,     # slightly more transparent
        smoothFactor = 0.5,    # reduces visual gaps
        highlightOptions = highlightOptions(
          weight = 2,
          color = "black",
          bringToFront = TRUE
        ),
        label = ~paste0(HYBAS_ID, ": ", round(pct_ghosts, 1), "% ghosts")
      ) %>%
      addLegend(
        "bottomright",
        pal = pal,
        values = ~pct_ghosts,
        title = "% molecular ghosts",
        labFormat = labelFormat(suffix = "%"),
        opacity = 0.7
      )
  })
  
  # Render DT table for clicked basin
  observeEvent(input$map_shape_click, {
    basin_id <- input$map_shape_click$id
    
    sp <- species_data %>%
      filter(HYBAS_ID == basin_id) %>%
      mutate(all = V5_12S + MiMamm_12S + Vences_16S) %>%
      arrange(all, sciname) %>%
      select(sciname, V5_12S, MiMamm_12S, Vences_16S)
    
    output$basin_table <- renderDT({
      datatable(
        sp,
        options = list(
          pageLength = 10,
          scrollY = "400px",
          scrollCollapse = TRUE
        ),
        rownames = FALSE
      )
    })
  })
}

# -------------------------------
# Run the app
# -------------------------------
shinyApp(ui, server)
