library(tidyverse)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

asvs <- read_tsv("outputs/final_tables/planktospace_metaB_v1_asvs.tsv.gz")
context <- read_tsv(
  "outputs/final_tables/planktospace_metaB_v1_context.tsv.gz"
) |>
  mutate(
    expedition = case_when(
      expedition == "europa" ~ "Tara Europa",
      expedition == "pacific" ~ "Tara Pacific"
    )
  )
counts <- read_tsv("outputs/final_tables/planktospace_metaB_v1_counts.tsv.gz")

# map

world <- ne_countries(scale = "medium", returnclass = "sf")

crsrobin <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

crsrobin_pacific <- "+proj=robin +lon_0=-180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

world_robinson_pacific <- st_break_antimeridian(world, lon_0 = 180) %>%
  st_transform(crs = crsrobin_pacific)

sites <- context %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    agr = "constant"
  )

ggplot(data = world) +
  geom_sf(colour = "black", fill = "gray85") +
  geom_sf(
    data = sites,
    aes(colour = expedition),
    size = 2,
    shape = 21,
    stroke = 1.2
  ) +
  coord_sf(crs = crsrobin) +
  theme_classic() +
  theme(
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 13),
    panel.grid.major = element_line(
      color = gray(.5),
      linetype = "dashed",
      linewidth = 0.5
    )
  )

ggsave("assets/img/sampling_map.png", width = 8, height = 6, dpi = 300)

# size fraction distribution

ggplot(context) +
  aes(
    x = as.factor(size_fraction_lower_threshold),
    fill = as.factor(size_fraction_upper_threshold)
  ) +
  geom_bar() +
  facet_wrap(~expedition) +
  scale_x_discrete("Size fraction lower threshold in µm") +
  scale_y_continuous("# samples") +
  scale_fill_discrete("Size fraction upper threshold in µm") +
  theme(legend.position = "bottom")

ggsave(
  "assets/img/size_fraction_distribution.png",
  width = 6,
  height = 4,
  dpi = 300
)


# taxonomic composition
asvs <-
  asvs |>
  mutate(
    domain = str_split_i(taxonomy, ";", 1),
    domain = case_when(
      domain %in% c("Archaea", "Bacteria", "Eukaryota") ~ domain,
      domain == "Eukaryota:plas" ~ "Chloroplast",
      is.na(domain) ~ "unassigned",
      .default = "others"
    ),
    domain = factor(
      domain,
      levels = c(
        "Archaea",
        "Bacteria",
        "Chloroplast",
        "Eukaryota",
        "others",
        "unassigned"
      )
    )
  )

taxo_comp <-
  asvs |>
  select(asv_id, domain) |>
  inner_join(counts, by = "asv_id") |>
  inner_join(
    select(context, sample, size_fraction_lower_threshold, expedition),
    by = "sample"
  ) |>
  group_by(domain, size_fraction_lower_threshold, expedition) |>
  summarise(nreads = sum(nreads), nasvs = n_distinct(asv_id))

a <- ggplot(taxo_comp) +
  aes(x = as.factor(size_fraction_lower_threshold), y = nreads, fill = domain) +
  geom_col() +
  facet_wrap(~expedition) +
  scale_x_discrete("Size fraction lower threshold in µm") +
  scale_y_continuous(
    "# reads",
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  scale_fill_brewer("Taxonomic group (domain)", palette = "Set1")

b <- ggplot(taxo_comp) +
  aes(x = as.factor(size_fraction_lower_threshold), y = nasvs, fill = domain) +
  geom_col() +
  facet_wrap(~expedition) +
  scale_x_discrete("Size fraction lower threshold in µm") +
  scale_y_continuous(
    "# ASVs",
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  scale_fill_brewer("Taxonomic group (domain)", palette = "Set1")


a +
  b +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "assets/img/domain_composition.png",
  width = 6,
  height = 6,
  dpi = 300
)
