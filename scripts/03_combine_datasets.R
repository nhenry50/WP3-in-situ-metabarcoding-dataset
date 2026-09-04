library(tidyverse)

### context

context <-
  bind_rows(
    read_tsv("outputs/subseted_tables/tara_europa_context_subseted.tsv.gz") |>
      mutate(expedition = "europa"),
    read_tsv("outputs/subseted_tables/tara_pacific_context_subseted.tsv.gz") |>
      mutate(expedition = "pacific")
  )

### counts

counts <-
  bind_rows(
    read_tsv("outputs/subseted_tables/tara_europa_counts_subseted.tsv.gz"),
    read_tsv("outputs/subseted_tables/tara_pacific_counts_subseted.tsv.gz")
  ) |>
  group_by(asv_id) |>
  filter(sum(nreads) >= 3, n_distinct(sample) >= 2) |>
  ungroup()

asvs_updated_stats <-
  counts |>
  group_by(asv_id) |>
  summarise(total = sum(nreads), spread = n_distinct(sample))

### asvs

asvs <-
  bind_rows(
    read_tsv("outputs/subseted_tables/tara_europa_asvs_subseted.tsv.gz"),
    read_tsv("outputs/subseted_tables/tara_pacific_asvs_subseted.tsv.gz")
  ) |>
  mutate(total = NULL, spread = NULL)

### in case taxonomic assignment differs, the highest confidence score is kept

asvs <-
  asvs |>
  mutate(confidence = ifelse(is.na(confidence), "0", confidence)) |>
  rowwise() |>
  mutate(
    sum_conf = str_split_1(confidence, ";") |>
      as.numeric() |>
      sum()
  ) |>
  group_by(asv_id) |>
  slice_max(order_by = sum_conf, with_ties = FALSE) |>
  ungroup() |>
  mutate(sum_conf = NULL)

asvs <-
  inner_join(asvs_updated_stats, asvs, by = "asv_id") |>
  arrange(desc(total))

# export

write_tsv(asvs, "outputs/final_tables/planktospace_metaB_v1_asvs.tsv.gz")
write_tsv(counts, "outputs/final_tables/planktospace_metaB_v1_counts.tsv.gz")
write_tsv(context, "outputs/final_tables/planktospace_metaB_v1_context.tsv.gz")
