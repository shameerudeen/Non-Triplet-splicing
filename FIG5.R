res_trip <- rnaseq_alt %>% filter(NMD == "Resistant",Group == "Triplets", sc_alt != sc_ref)
res_ntp <- rnaseq_alt %>% filter(NMD == "Resistant",Group == "Non-Triplets")
resis <- bind_rows(res_ntp,res_trip)

resis <- resis %>% mutate(mechanism = case_when( distance_sc_alt < 80 ~ "EJC",TRUE ~ "Unknown"))
resis %>% group_by(mechanism) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100, 0))
1 EJC         112   210      53
2 Unknown      98   210      47
resis %>% group_by(exon_after_sc_alt) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100, 0))
exon_after_sc_alt count total percent

1 0                    76   210      36
2 1                    53   210      25
3 2                    19   210       9
4 3                    12   210       6
5 4                    17   210       8
6 5                     4   210       2
7 6                     9   210       4
8 7                     5   210       2
9 8                     3   210       1
10 9                     1   210       0
11 10                    6   210       3
12 12                    3   210       1
13 13                    2   210       1


# Step 1: Summarize count of events
resis_summary <- resis %>%
  group_by(exon_after_sc_alt) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    exon_after_sc_alt = factor(exon_after_sc_alt, levels = sort(unique(exon_after_sc_alt)))
  )

# Step 2: Plot count data
ggplot(resis_summary, aes(x = exon_after_sc_alt, y = count)) +
  geom_bar(stat = "identity", fill = "#00BFC4", color = "black", width = 0.7) +
  labs(
    x = "Number of exon junctions after the PTC",
    y = "Number of events"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(face = "bold", size = 11)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))



dual_code <- resis %>% filter(Group == "Non-Triplets",dual_coding > 60)

c_ter <- resis %>% filter(!(Group == "Non-Triplets" & dual_coding > 60))
c_ter %>% group_by(category) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100, 0))

