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

sample_provenance <-
  read_tsv(
    "data/external_tables/tara_pacific_context/TARA-PACIFIC_samples-provenance_20220131d-1.tsv.gz",
    skip = 1
  )

sample_provenance <-
  sample_provenance |>
  filter(
    `sampling-environment_feature_label` %in%
      c("O-SRF", "I-SRFa", "I-SRFb", "S-SRF")
  ) |>
  select(
    sample = `sample-id_source`,
    id_biosamples = `sample-id_biosamples`,
    size_fraction_lower_threshold = `sample-material_size-fraction_lower-threshold_micrometre`,
    size_fraction_upper_threshold = `sample-material_size-fraction_upper-threshold_micrometre`,
    datetime_utc_start = `sampling-event_datetime-utc_start_yyyy-mm-ddThh:mm:ssZ00`,
    datetime_utc_end = `sampling-event_datetime-utc_end_yyyy-mm-ddThh:mm:ssZ00`,
    latitude_start = `sampling-event_latitude_start_dd.dddddd`,
    latitude_end = `sampling-event_latitude_end_dd.dddddd`,
    longitude_start = `sampling-event_longitude_start_ddd.dddddd`,
    longitude_end = `sampling-event_longitude_end_ddd.dddddd`,
    depth = `sampling-event_depth-below-sea-surface_nominal`
  ) |>
  mutate(
    sample = str_remove(sample, "TARA_") |> str_replace("-", "_"),
    size_fraction_lower_threshold = replace(
      size_fraction_lower_threshold,
      size_fraction_lower_threshold == "<",
      NA
    ),
    size_fraction_upper_threshold = replace(
      size_fraction_upper_threshold,
      size_fraction_upper_threshold == ">",
      NA
    ),
    depth = str_remove(depth, " m")
  ) |>
  rowwise() |>
  mutate(
    datetime_utc = mean(c(datetime_utc_start, datetime_utc_end)),
    latitude = mean(c(latitude_start, latitude_end)),
    longitude = mean(c(longitude_start, longitude_end))
  ) |>
  ungroup() |>
  select(-c(ends_with("_start"), ends_with("_end")))

#####################################################################
# import and subset count data
#####################################################################

counts <-
  read_tsv(
    "data/external_tables/tara_pacific_asvs/tara_pacific_ssu_v4v5_dada2_counts.tsv.gz"
  )

counts <-
  sample_provenance |>
  select(sample, id_biosamples) |>
  distinct() |>
  inner_join(counts, by = "sample") |>
  mutate(sample = NULL) |>
  rename(sample = id_biosamples)

# subset sample context based on available counts

sample_provenance <-
  sample_provenance |>
  filter(id_biosamples %in% counts$sample) |>
  mutate(sample = NULL) |>
  rename(sample = id_biosamples)


#####################################################################
# import and subset asv taxonomy
#####################################################################

asvs <- read_tsv(
  "data/external_tables/tara_pacific_asvs/tara_pacific_ssu_v4v5_dada2_asvs.tsv.gz"
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

write_tsv(asvs, "outputs/subseted_tables/tara_pacific_asvs_subseted.tsv.gz")
write_tsv(counts, "outputs/subseted_tables/tara_pacific_counts_subseted.tsv.gz")
write_tsv(
  sample_provenance,
  "outputs/subseted_tables/tara_pacific_context_subseted.tsv.gz"
)
