# Spatial distribution maps for Phase 3 submission
# Full epidemic (560 cases, 22 Dec 2025 – 7 Apr 2026)
# Output: submission/spatial_maps_phase3.pdf

library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
here::i_am("hpai-challenge/analysis/phase-3-spatial-maps.R")

# Load data
cases <- read_csv(here::here("hpai-challenge/data/phase-3/cases.csv"),
                   show_col_types = FALSE) |>
  mutate(date_confirmed = ymd(date_confirmed))

population <- read_csv(here::here("hpai-challenge/data/canonical/population.csv"),
                        show_col_types = FALSE)

cases_full <- cases |>
  left_join(population, by = "farm_id") |>
  mutate(
    epidemic_day = as.numeric(date_confirmed - ymd("2025-12-22")),
    phase = case_when(
      date_confirmed <= ymd("2026-01-13") ~ "Phase 1",
      date_confirmed <= ymd("2026-02-13") ~ "Phase 2",
      TRUE ~ "Phase 3"
    )
  )

theme_set(
  theme_minimal(base_size = 11) +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
)

# Map 1: Outbreaks by species
p1 <- ggplot() +
  geom_point(data = population, aes(x = x, y = y),
             colour = "grey80", size = 0.3, alpha = 0.5) +
  geom_point(data = cases_full,
             aes(x = x, y = y, colour = species, shape = species),
             size = 2, alpha = 0.8) +
  scale_colour_manual(
    values = c("duck" = "#1f77b4", "chicken" = "#ff7f0e"),
    labels = str_to_title, name = "Species"
  ) +
  scale_shape_manual(
    values = c("duck" = 16, "chicken" = 17),
    labels = str_to_title, name = "Species"
  ) +
  coord_equal() +
  labs(title = paste0("A) Outbreaks by species (n = ", nrow(cases_full), ")")) +
  theme(legend.position = "bottom")

# Map 2: Temporal progression
p2 <- ggplot() +
  geom_point(data = population, aes(x = x, y = y),
             colour = "grey90", size = 0.2, alpha = 0.4) +
  geom_point(data = cases_full |> arrange(epidemic_day),
             aes(x = x, y = y, colour = epidemic_day),
             size = 2, alpha = 0.85) +
  scale_colour_viridis_c(
    name = "Days since\nfirst case", option = "inferno", direction = -1
  ) +
  coord_equal() +
  labs(title = "B) Temporal progression") +
  theme(legend.position = "bottom")

# Map 3: By production type
p3 <- ggplot() +
  geom_point(data = population, aes(x = x, y = y),
             colour = "grey90", size = 0.2, alpha = 0.4) +
  geom_point(data = cases_full,
             aes(x = x, y = y, colour = production),
             size = 2, alpha = 0.8) +
  scale_colour_brewer(
    palette = "Set1",
    labels = \(x) str_replace_all(x, "_", " ") |> str_to_title(),
    name = "Production type"
  ) +
  coord_equal() +
  labs(title = "C) Outbreaks by production type") +
  theme(legend.position = "bottom")

# Save
pdf(here::here("hpai-challenge/submission/spatial_maps_phase3.pdf"),
    width = 12, height = 6)
print(p1 + p2 + p3 + plot_layout(ncol = 3))
dev.off()

cat("Saved: submission/spatial_maps_phase3.pdf\n")
