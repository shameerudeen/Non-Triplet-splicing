
# process cassette exons from JUM

cas_dtl <- fread("cas_detailed.txt",fill=TRUE)

cas_dtl <- cas_dtl %>%
  mutate(
    parts = strsplit(AS_event_ID, "_"),
    chr = sapply(parts, `[`, 1),
    strand = sapply(parts, `[`, 2),
    start = as.integer(sapply(parts, `[`, 4)) + 1,
    end = as.integer(sapply(parts, `[`, 5))
  ) %>%
  dplyr::select(-parts)


cas_dtl_clean <- cas_dtl %>%
  group_by(Gene, start, end, strand) %>%
  mutate(row_index = row_number()) %>%
  mutate(deltaPSI = ifelse(row_index == 3,
                           deltaPSI[row_index == 4],
                           deltaPSI)) %>%
  filter(row_index %in% c(2, 3)) %>%
  ungroup()

cas_dtl_clean <- cas_dtl_clean %>%
  mutate(across(c(percentage_usage.control1,
                  percentage_usage.control2,
                  percentage_usage.control3,
                  percentage_usage.upf1,
                  percentage_usage.upf2,
                  percentage_usage.upf3),
                ~ as.numeric(gsub("%", "", .x))))

library(dplyr)

cas_dtl_final <- cas_dtl_clean %>%
  group_by(Gene, start, end, strand) %>%
  mutate(
    avgwtPSI = mean(c(
      percentage_usage.control1,
      percentage_usage.control2,
      percentage_usage.control3
    )),
    avgUPF1PSI = mean(c(
      percentage_usage.upf1,
      percentage_usage.upf2,
      percentage_usage.upf3
    ))
  ) %>%
  mutate(
    avgwtPSI = ifelse(row_number() == 2,
                      mean(c_across(starts_with("percentage_usage.control"))),
                      NA),
    avgUPF1PSI = ifelse(row_number() == 2,
                        mean(c_across(starts_with("percentage_usage.upf"))),
                        NA)
  ) %>%
  filter(row_number() == 2) %>%
  ungroup()

library(dplyr)

cassette <- cas_dtl_final %>%
  filter(
    !is.na(deltaPSI), !is.infinite(deltaPSI),
    !is.na(avgwtPSI), !is.infinite(avgwtPSI),
    !is.na(avgUPF1PSI), !is.infinite(avgUPF1PSI)
  ) %>%
  select(-row_index)

cassette <- cassette %>% mutate(deltaPSI = -deltaPSI)
cassette <- cassette %>% select(Gene,AS_event_ID,`BH_adjusted_p-values`,chr,strand,start,end,avgwtPSI,avgUPF1PSI,deltaPSI) %>% 
  mutate(category= "cassette")

cassette <- cassette %>% rename (Significance = `BH_adjusted_p-values`)



# A3SS

a3_dtl <- fread("A3SS_detailed_with_gene.txt",fill=TRUE)


a3_dtl <- a3_dtl %>%
  mutate(
    parts = strsplit(AS_event_ID, "_"),
    chr = sapply(parts, `[`, 1),
    strand = sapply(parts, `[`, 2),
    start = case_when(
      strand == "-" ~ as.integer(sapply(parts, `[`, 3)) + 1,
      strand == "+" ~ as.integer(sapply(parts, `[`, 4)) + 2
    ),
    end = case_when(
      strand == "-" ~ as.integer(sapply(parts, `[`, 4)),
      strand == "+" ~ as.integer(sapply(parts, `[`, 5)) + 1
    )
  ) %>%
  dplyr::select(-parts)

a3_dtl <- a3_dtl %>%
  filter(lengths(strsplit(AS_event_ID, "_")) == 5)



a3_dtl_minus <- a3_dtl %>%
  filter(strand == "-") %>%
  group_by(Gene, chr, strand, start, end) %>%
  slice_min(sub_junction_start_coor, n = 1, with_ties = FALSE) %>%
  ungroup()

a3_dtl_plus <- a3_dtl %>%
  filter(strand == "+") %>%
  group_by(Gene, chr, strand, start, end) %>%
  slice_min(sub_junction_end_coor, n = 1, with_ties = FALSE) %>%
  ungroup()

a3ss_dtl <- bind_rows(a3_dtl_minus, a3_dtl_plus) 



a3ss_dtl <- a3ss_dtl %>%
  mutate(
    across(
      c(percentage_usage.control1,
        percentage_usage.control2,
        percentage_usage.control3,
        percentage_usage.upf1,
        percentage_usage.upf2,
        percentage_usage.upf3),
      ~ as.numeric(gsub("%", "", .x))
    ),
    avgwtPSI   = rowMeans(across(c(percentage_usage.control1,
                                   percentage_usage.control2,
                                   percentage_usage.control3)),
                          na.rm = TRUE),
    avgUPF1PSI = rowMeans(across(c(percentage_usage.upf1,
                                   percentage_usage.upf2,
                                   percentage_usage.upf3)),
                          na.rm = TRUE)
  )

a3ss <- a3ss_dtl %>%
  filter(
    !is.na(deltaPSI), !is.infinite(deltaPSI),
    !is.na(avgwtPSI), !is.infinite(avgwtPSI),
    !is.na(avgUPF1PSI), !is.infinite(avgUPF1PSI)
  ) 

a3ss <- a3ss %>% dplyr::select(Gene,AS_event_ID,`BH_adjusted_p-values`,chr,strand,start,end,avgwtPSI,avgUPF1PSI,deltaPSI) %>% 
  mutate(category= "a3ss")

a3ss <- a3ss %>% dplyr::rename (Significance = `BH_adjusted_p-values`)



# A5SS


a5_dtl <- fread("A5SS_detailed.txt",fill=TRUE)

a5_dtl <- a5_dtl %>%
  mutate(
    parts = strsplit(AS_event_ID, "_"),
    chr = sapply(parts, `[`, 1),
    strand = sapply(parts, `[`, 2),
    start = case_when(
      strand == "-" ~ as.integer(sapply(parts, `[`, 4)) + 2,
      strand == "+" ~ as.integer(sapply(parts, `[`, 3)) + 1
    ),
    end = case_when(
      strand == "-" ~ as.integer(sapply(parts, `[`, 5)) + 1,
      strand == "+" ~ as.integer(sapply(parts, `[`, 4))
    )
  ) %>%
  dplyr::select(-parts)


a5_dtl <- a5_dtl %>%
  filter(lengths(strsplit(AS_event_ID, "_")) == 5)

a5_dtl_minus <- a5_dtl %>%
  filter(strand == "-") %>%
  group_by(Gene, chr, strand, start, end) %>%
  slice_min(sub_junction_end_coor, n = 1, with_ties = FALSE) %>%
  ungroup()

a5_dtl_plus <- a5_dtl %>%
  filter(strand == "+") %>%
  group_by(Gene, chr, strand, start, end) %>%
  slice_min(sub_junction_start_coor, n = 1, with_ties = FALSE) %>%
  ungroup()

a5ss_dtl <- bind_rows(a5_dtl_minus, a5_dtl_plus) 



a5ss_dtl <- a5ss_dtl %>%
  mutate(
    across(
      c(percentage_usage.control1,
        percentage_usage.control2,
        percentage_usage.control3,
        percentage_usage.upf1,
        percentage_usage.upf2,
        percentage_usage.upf3),
      ~ as.numeric(gsub("%", "", .x))
    ),
    avgwtPSI   = rowMeans(across(c(percentage_usage.control1,
                                   percentage_usage.control2,
                                   percentage_usage.control3)),
                          na.rm = TRUE),
    avgUPF1PSI = rowMeans(across(c(percentage_usage.upf1,
                                   percentage_usage.upf2,
                                   percentage_usage.upf3)),
                          na.rm = TRUE)
  )

a5ss <- a5ss_dtl %>%
  filter(
    !is.na(deltaPSI), !is.infinite(deltaPSI),
    !is.na(avgwtPSI), !is.infinite(avgwtPSI),
    !is.na(avgUPF1PSI), !is.infinite(avgUPF1PSI)
  ) 

a5ss <- a5ss %>% dplyr::select(Gene,AS_event_ID,`BH_adjusted_p-values`,chr,strand,start,end,avgwtPSI,avgUPF1PSI,deltaPSI) %>% 
  mutate(category= "a5ss")

a5ss <- a5ss %>% dplyr::rename (Significance = `BH_adjusted_p-values`)





