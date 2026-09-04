library(tidyverse)

#####################################################################
# function to make a consensual taxonomy between silva and pr2,
# with priority of assignation to pr2 for eukaryotes and silva
# for prokaryotes
#####################################################################

make_consensus_taxo <- function(taxonomy) {
  taxonomy |>
    mutate(
      consensus_dada2_taxonomy = case_when(
        str_detect(pr2_dada2_taxonomy, "^Eukaryota") ~ pr2_dada2_taxonomy,
        str_detect(silva_dada2_taxonomy, "Mitochondria") ~ "Eukaryota:mito",
        str_detect(silva_dada2_taxonomy, "^Eukaryota|Chloroplast") ~ "conflict",
        str_detect(
          silva_dada2_taxonomy,
          "^Bacteria|^Archaea"
        ) ~ silva_dada2_taxonomy,
        pr2_dada2_taxonomy == "" | silva_dada2_taxonomy == "" ~ "unassigned",
        .default = NA
      ),
      consensus_dada2_confidence = case_when(
        str_detect(pr2_dada2_taxonomy, "^Eukaryota") ~ pr2_dada2_confidence,
        str_detect(
          silva_dada2_taxonomy,
          "Mitochondria"
        ) ~ silva_dada2_confidence,
        str_detect(silva_dada2_taxonomy, "^Eukaryota|Chloroplast") ~ NA,
        str_detect(
          silva_dada2_taxonomy,
          "^Bacteria|^Archaea"
        ) ~ silva_dada2_confidence,
        pr2_dada2_taxonomy == "" | silva_dada2_taxonomy == "" ~ NA,
        .default = NA
      )
    )
}

#####################################################################
# import and subset contextual data
#####################################################################

context <- read_tsv(
  "data/external_tables/tara_europa_asvs/trec-soil_sediment_waters-ssuv4v5_biosamples.tsv.gz"
)

str(context)

context |>
  pull(sample_alias)

context <-
  context |>
  filter(str_detect(sample_alias, "_water_")) |>
  select(
    sample = accession,
    size_fraction_lower_threshold = `size-fraction lower threshold`,
    size_fraction_upper_threshold = `size-fraction upper threshold`,
    datetime_utc = `collection date`,
    latitude = `geographic location (latitude)`,
    longitude = `geographic location (longitude)`,
    depth
  ) |>
  mutate(
    sample = str_remove(sample, "TARA_") |> str_replace("-", "_"),
    size_fraction_lower_threshold = replace(
      size_fraction_lower_threshold,
      size_fraction_lower_threshold == "not applicable",
      NA
    ),
    size_fraction_upper_threshold = replace(
      size_fraction_upper_threshold,
      size_fraction_upper_threshold == "not applicable",
      NA
    )
  )

#####################################################################
# import and subset count data
#####################################################################

counts <-
  read_tsv(
    "data/external_tables/tara_europa_asvs/trec-soil_sediment_waters-ssuv4v5_dada2_counts.tsv.gz"
  ) |>
  filter(sample %in% context$sample)

#####################################################################
# import and subset asv taxonomy
#####################################################################

asvs <- read_tsv(
  "data/external_tables/tara_europa_asvs/trec-soil_sediment_waters-ssuv4v5_dada2_asvs.tsv.gz"
) |>
  filter(asv_id %in% counts$asv_id)

asvs_updated_stats <-
  counts |>
  group_by(asv_id) |>
  summarise(total = sum(nreads), spread = n_distinct(sample))

asvs <-
  asvs |>
  make_consensus_taxo() |>
  mutate(total = NULL, spread = NULL) |>
  inner_join(asvs_updated_stats, by = "asv_id") |>
  select(
    asv_id,
    total,
    spread,
    taxonomy = consensus_dada2_taxonomy,
    confidence = consensus_dada2_confidence,
    sequence
  ) |>
  arrange(desc(total))

#####################################################################
# export subseted tables
#####################################################################

write_tsv(asvs, "outputs/subseted_tables/tara_europa_asvs_subseted.tsv.gz")
write_tsv(counts, "outputs/subseted_tables/tara_europa_counts_subseted.tsv.gz")
write_tsv(
  context,
  "outputs/subseted_tables/tara_europa_context_subseted.tsv.gz"
)
