cassette_clean$sc_alt <- as.integer(cassette_clean$sc_alt)
a3_clean$sc_alt <- as.integer(a3_clean$sc_alt)
a5_clean$sc_alt <- as.integer(a5_clean$sc_alt)

overall_clean <- bind_rows(cassette_clean, a3_clean, a5_clean)
overall_clean<-overall_clean %>% left_join(EXP %>% dplyr::select(Gene, log2FC_correct, padj), by = "Gene")
overall_clean <- overall_clean %>% left_join(smg_2 %>% dplyr::select(Gene, AS_event_ID, smg_2_delPSI, smg_2_wtPSI, smg_2PSI),by = c("Gene", "AS_event_ID"))
write.xlsx(overall_clean, file = "overall_clean.xlsx", overwrite = TRUE)
rnaseq_alt <- overall_clean %>% filter((avgwtPSI >= 0.10 & avgwtPSI <= 0.90) | (Abs_deltaPSI >= 0.10 & Significance < 0.05))
wt <- overall_clean %>% filter(avgwtPSI >= 0.10, avgwtPSI <= 0.90)
wt %>% group_by(Group) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))

#1507+ 7= 1514 events
smg_alt <- overall_clean %>% filter(!(avgwtPSI >= 0.10 & avgwtPSI <= 0.90) & Abs_deltaPSI >= 0.10 &  Significance < 0.05)
smg_alt %>% group_by(category) %>% summarise(n = n())

write.xlsx(rnaseq_alt, file = "rnaseq_alt.xlsx", overwrite = TRUE)


rnaseq_alt %>% group_by(category) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 a3ss       975  1885    51.7
2 a5ss       468  1885    24.8
3 cassette   442  1885    23.4

library(ggplot2)
library(dplyr)
library(ggrepel)
library(svglite)


# Data
event_counts <- data.frame(Event = c("Cassette", "A3SS", "A5SS", "MXE"),Count = c(442, 975, 468, 7))

# Desired order
event_order <- c("Cassette", "A3SS", "A5SS", "MXE")

event_counts <- event_counts %>% mutate(fraction = Count / sum(Count),ymax = cumsum(fraction),ymin = c(0, head(ymax, -1)),label_position = (ymax + ymin) / 2,
    label_text = Event,legend_text = paste0(Event, " (", Count, ")"))

# Explicit factor ordering for both Event and legend_text
event_counts$Event <- factor(event_counts$Event, levels = event_order)
event_counts$legend_text <- factor(event_counts$legend_text,
                                   levels = paste0(event_order, " (", event_counts$Count, ")"))

# Define colors following the same order
event_colors <- c("Cassette" = "#F8766D", "A3SS" = "darkseagreen4","A5SS" = "deeppink4","MXE" = "lemonchiffon3")

legend_colors <- setNames(event_colors, levels(event_counts$legend_text))

# Plot
A <- ggplot(event_counts, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0)) +
  geom_rect(aes(fill = legend_text), color = "black") +
  coord_polar(theta = "y") +
  xlim(c(-0.5, 1.5)) +
  scale_fill_manual(values = legend_colors, drop = FALSE) +
  geom_text( aes(x = 0.5, y = label_position, label = label_text),size = 6, fontface = "bold", color = "black"
  ) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(face = "bold", size = 20),
    legend.key.height = unit(1.5, "cm"),   
    legend.key.width  = unit(1.2, "cm"))


ggsave("1A.svg",A , width = 7, height = 7, units = "in", device = svglite::svglite)



rnaseq_alt %>% group_by(Group) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))
rnaseq_alt %>% group_by(category) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))

# Non-Triplets   684/  1892      36
# Triplets      1208/ 1892    64

rnaseq_alt %>% group_by(category,Group) %>% summarise(count = n(), .groups = "drop") %>% group_by(category) %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 a3ss     Non-Triplets   278   975      29
2 a3ss     Triplets       697   975      71
3 a5ss     Non-Triplets   272   468      58
4 a5ss     Triplets       196   468      42
5 cassette Non-Triplets   134   442      30
6 cassette Triplets       308   442      70

Per_tp <- data.frame(
  DataFrame = c("Overall", "Cassette", "A3SS", "A5SS","MXE"),
  Triplet_Percentage = c(64,70,71,42,100)
)

# Set the order of the DataFrame column
Per_tp$DataFrame <- factor(Per_tp$DataFrame, levels = c("Overall", "Cassette", "A3SS", "A5SS","MXE"))

# Percent Triplet
B<- ggplot(Per_tp, aes(x = DataFrame, y = Triplet_Percentage, fill = DataFrame)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +  # Increase width to reduce gaps
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "darkseagreen4","deeppink4", "lemonchiffon3")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Add panel border
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 24),  # Rotate x-axis labels
    axis.text.y = element_text(face = "bold", size = 20),  # Bold y-axis text
    axis.ticks = element_line(color = "black", size = 1),  # X-axis ticks
    axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 10)),  # Y-axis title spacing
    axis.title.x = element_text(face = "bold", size = 24),  # Bold x-axis title   
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # No title
  ) +
  labs(x = "",  y = "Percent Triplet" ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Bars start from x-axis
  scale_x_discrete(expand = c(0.1, 0.2))  # Reduce space between bars

ggsave("1B.svg", B, width = 4.5, height = 6, units = "in", device = svglite::svglite)

#FRAMR DISRUPTION
ntp <- rnaseq_alt %>% filter(Group== "Non-Triplets")

ntp %>% group_by(NMD) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 Resistant   190   684      28
2 Sensitive   494   684      72

ntp %>% group_by(category,NMD) %>% summarise(count = n(), .groups = "drop") %>% group_by(category) %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 a3ss     Resistant    82   278      29
2 a3ss     Sensitive   196   278      71
3 a5ss     Resistant    52   272      19
4 a5ss     Sensitive   220   272      81
5 cassette Resistant    56   134      42
6 cassette Sensitive    78   134      58

trp <- rnaseq_alt %>% filter(Group== "Triplets")
trp %>% group_by(NMD) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 Resistant   978  1201      81
2 Sensitive   223  1201      19
trp %>% group_by(category,NMD) %>% summarise(count = n(), .groups = "drop") %>% group_by(category) %>% mutate(total = sum(count),percent = round((count / total) * 100))

nmd <- data.frame(
  DataFrame = c("Overall", "Cassette", "A3SS", "A5SS"),Non_Triplet_Percentage = c(72,58,71,81))

nmd$DataFrame <- factor(nmd$DataFrame, levels = c("Overall", "Cassette", "A3SS", "A5SS"))

# Plot
C <- ggplot(nmd, aes(x = DataFrame, y = Non_Triplet_Percentage, fill = DataFrame)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D", "darkseagreen4","deeppink4")) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
        panel.border = element_rect(color = "black", fill = NA),  # Add panel border
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 24),  # Rotate x-axis labels
        axis.text.y = element_text(face = "bold", size = 20),  # Bold y-axis text
        axis.ticks = element_line(color = "black", size = 1),  # X-axis ticks
        axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 10)),  # Y-axis title spacing
        axis.title.x = element_text(face = "bold", size = 24),  # Bold x-axis title   
        legend.position = "none",  # Remove legend
        plot.title = element_blank()  # No title
  ) +
  labs(x = "",  y = "Percent NMD Sensitive" ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +  
  scale_x_discrete(expand = c(0.1, 0.2))  # Reduce space between bars

ggsave("1C.svg", C, width = 4, height = 6, units = "in", device = svglite::svglite)

# splice site selection
calculate_percent_nontriplet_plot_and_table <- function(data) {
  library(dplyr)
  library(ggplot2)
  
  # Define PSI bins
  psi_groups <- c("1-20%", "20-40%", "40-60%", "60-80%", "80-100%")
  psi_breaks <- c(0, 0.20, 0.40, 0.60, 0.80, 1.00)
  
  # Function to calculate percent Non-Triplets per PSI bin
  calculate_percent <- function(df, psi_col) {
    df %>%
      mutate(PSI_Bin = cut(
        !!sym(psi_col),
        breaks = psi_breaks,
        labels = psi_groups,
        include.lowest = TRUE
      )) %>%
      group_by(PSI_Bin) %>%
      summarise(
        total = n(),
        nontriplet_count = sum(Group == "Non-Triplets", na.rm = TRUE),
        PercentNonTriplet = round((nontriplet_count / total) * 100, 1),
        .groups = "drop"
      )
  }
  
  # Calculate for WT and smg-1
  wt_data   <- calculate_percent(data, "avgwtPSI")  %>% mutate(Type = "WT")
  smg1_data <- calculate_percent(data, "avgsmgPSI") %>% mutate(Type = "smg-1")
  
  final_data <- bind_rows(wt_data, smg1_data) %>%
    mutate(Type = factor(Type, levels = c("WT", "smg-1")))  # enforce order
  
  # Plot
  plot <- ggplot(final_data, aes(x = PSI_Bin, y = PercentNonTriplet, color = Type, group = Type)) +
    geom_point(size = 4) +
    geom_line(linewidth = 1.2) +
    labs(
      x = "Percent Spliced In",
      y = "Percent Non Triplet",
      color = NULL
    ) +
    ylim(0, 70) +
    scale_color_manual(
      values = c("WT" = "#F8766D", "smg-1" = "#00BFC4"),
      labels = c(
        "WT"    = expression(bold("Wild Type")),
        "smg-1" = expression(bolditalic("smg-1"))
      )
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 20),
      axis.text.y = element_text(face = "bold", size = 20),
      axis.ticks = element_line(color = "black", size = 1),
      axis.title.y = element_text(face = "bold", size = 24, margin = margin(r = 10)),
      axis.title.x = element_text(face = "bold", size = 24),
      legend.text = element_text(face = "bold", size = 24),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width  = unit(1.2, "cm")
    ) +
    labs(title = "")
  
  list(
    plot = plot,
    table = final_data %>%
      dplyr::select(Type, PSI_Bin, PercentNonTriplet, total, nontriplet_count)
  )
}



# Run the function
result <- calculate_percent_nontriplet_plot_and_table(overall_clean)

# View the plot
D <- result$plot

ggsave("1D.svg", D, width = 7, height = 6, units = "in", device = svglite::svglite)

# View the table
print(result$table)

1 WT    1-20%                62.5   858              536
2 WT    20-40%               17.5   354               62
3 WT    40-60%               14.2   295               42
4 WT    60-80%               19.8   373               74
5 WT    80-100%              65.1   870              566
6 smg-1 1-20%                54.4   759              413
7 smg-1 20-40%               36.8   454              167
8 smg-1 40-60%               31.5   352              111
9 smg-1 60-80%               37.3   413              154
10 smg-1 80-100%              56.3   772              435




rnaseq_alt %>% group_by(NMD) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 Resistant  1168  1885      62
2 Sensitive   717  1885      38

sens <- rnaseq_alt %>% filter(NMD == "Sensitive", sc_alt != sc_ref)
sens %>% group_by(Group) %>% summarise(count = n(), .groups = "drop") %>% mutate(total = sum(count),percent = round((count / total) * 100))
1 Non-Triplets   494   546      90
2 Triplets        52   546      10










library(dplyr)

volcano <- rnaseq_alt %>%
  filter(Group == "Non-Triplets") %>%
  mutate(deltaPSI_percent = round(deltaPSI * 100, 0)) %>%
  dplyr::select(-any_of("NMD")) %>%
  mutate(NMD = ifelse(abs(deltaPSI_percent) > 10, "Sensitive", "Resistant"))

library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)   

# Define genes of interest
genes_of_interest <- c("dhhc-4", "fubl-3", "tnt-3", "fubl-1")

E <- ggplot(volcano, aes(x = deltaPSI_percent, y = -log10(Significance))) +
  geom_vline(xintercept = c(-10, 10), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_jitter(data = subset(volcano, !(Gene %in% genes_of_interest)),
              aes(color = NMD),
              width = 0.5, height = 0.5, size = 2, alpha = 1) +
  geom_jitter(data = subset(volcano, Gene %in% genes_of_interest),
              color = "black", width = 0.5, height = 0.5, size = 2, alpha = 1,show.legend = FALSE) +
  geom_text_repel(data = subset(volcano, Gene %in% genes_of_interest), aes(label = Gene),
                  size = 9, family = "Helvetica Neue", fontface = "bold.italic",
                  max.overlaps = 20, box.padding = 0.5, point.padding = 0.3) +
  scale_color_manual(name = "Splicing event",
                     values = c("Sensitive" = "#F8766D", "Resistant" = "#00BFC4"),
                     labels = c("Sensitive" = "NMD Sensitive", "Resistant" = "NMD Resistant")) +
  scale_y_continuous(name = expression(bold(-log[10]("Significance"))),
                     limits = c(0, 120), breaks = seq(0, 120, by = 30)
  ) +
  scale_x_continuous(
    name = expression(bold(Delta) * bold("Percent Spliced In (") * bolditalic("smg-1") * bold(" vs WT)")),
    limits = c(-90, 90), breaks = seq(-90, 90, by = 30)
  ) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 7)),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        legend.key.height = unit(1.5, "cm"),   legend.key.width  = unit(1.2, "cm"), 
        panel.border = element_rect(color = "black", fill = NA),panel.grid = element_blank())



ggsave("E.svg", E, width = 12, height = 6, units = "in", device = svglite::svglite)


library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# ---- Data ----
df1 <- data.frame(
  Event = c("Cassette", "Cassette","A3SS", "A3SS", "A5SS", "A5SS"),
  Group = c("Non-Triplets", "Triplets", "Non-Triplets", "Triplets", "Non-Triplets", "Triplets"),
  Count = c(78, 28,196, 9, 220, 15)
)

df2 <- data.frame(
  Event = c("Cassette", "Cassette","A3SS", "A3SS", "A5SS", "A5SS"),
  Group = c("Non-Triplets", "Triplets", "Non-Triplets", "Triplets", "Non-Triplets", "Triplets"),
  Count = c(38, 9,60, 6, 38, 5)
)

df3 <- data.frame(
  Event = c("Cassette", "Cassette","A3SS", "A3SS", "A5SS", "A5SS"),
  Group = c("Non-Triplets", "Triplets", "Non-Triplets", "Triplets", "Non-Triplets", "Triplets"),
  Count = c(18, NA,22, NA, 14, NA)   # Triplets = NA
)

# Factor order + colors
event_levels <- c("Cassette", "A3SS", "A5SS")
event_colors <- c("Cassette" = "#F8766D",
                  "A3SS" = "darkseagreen4",
                  "A5SS" = "deeppink4")

# ---- Function to make a horizontal bar plot ----
make_plot <- function(df, show_title = NULL, show_legend = FALSE) {
  df$Event <- factor(df$Event, levels = event_levels)  # keep legend order
  
  ggplot(df, aes(x = Group, y = Count, fill = Event)) +
    geom_bar(stat = "identity", width = 0.5, na.rm = TRUE,
             position = position_stack(reverse = TRUE)) +  # reverse stack only
    scale_fill_manual(values = event_colors, breaks = event_levels) +
    labs(x = NULL, y = NULL, title = show_title) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  
      panel.border = element_rect(color = "black", fill = NA),  
      axis.text.x = element_text(face = "bold", size = 20),  
      axis.text.y = element_text(face = "bold", size = 24),  
      axis.ticks = element_line(color = "black", size = 1),  
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      legend.title = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(face = "bold", size = 18),
      legend.key.height = unit(0.7, "cm"),   
      legend.key.width  = unit(0.8, "cm")
    ) +
    coord_flip()
}

# ---- Make plots ----
F1 <- make_plot(df1, show_title = "Splicing events", show_legend = FALSE)
F2 <- make_plot(df2, show_title = NULL, show_legend = FALSE)
F3 <- make_plot(df3, show_title = NULL, show_legend = TRUE)

# ---- Combine vertically ----
combined <- gridExtra::grid.arrange(F1, F2, F3, ncol = 1)

# ---- Save ----
ggsave("F.svg", combined, width = 6, height = 7, units = "in", device = svglite::svglite)

# Violin plot for CASSETTE dataset with NMD-colored points
ggplot(CASSETTE, aes(x=Group, y=average_phyloP, fill=Group)) +
  geom_violin(trim=FALSE, alpha=0.7) +  # Violin plot for density
  geom_point(aes(color=NMD), position=position_jitter(width=0.2), size=1, alpha=0.7) +  # Jittered points
  geom_hline(yintercept = 1, linetype = "dotted", color = "blue", size = 1) +  # Add dotted blue line at y=1
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add panel border
    axis.text.x = element_text(hjust = 0.5, face = "bold", size = 12),  # X-axis labels rotated and bold
    axis.text.y = element_text(face = "bold", size = 10),  # Y-axis text bold
    axis.ticks.x = element_line(color = "black", size = 0.5),  # Add x-axis ticks
    axis.ticks.y = element_line(color = "black", size = 0.5),  # Add y-axis ticks
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 10)),  # Adjust y-axis title
    axis.title.x = element_blank(),  # Remove x-axis title
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # No title
  ) +
  labs(
    x = "",  # X-axis label removed
    y = "Conservation score"  # Updated y-axis label
  ) +
  scale_fill_manual(values = c("Triplets" = "cadetblue3", "Non-Triplets" = "cornsilk3")) +  # Custom violin colors
  scale_color_manual(values = c("Sensitive"="red", "Resistant"="black"))  # Color points by NMD

library(dplyr)

# Filter CASSETTE for Non-Triplets, Sensitive, and average_phyloP > 2
consNTP_CASS <-CASSETTE%>%
  filter(Group == "Non-Triplets", NMD == "Sensitive", average_phyloP < 1)  %>%
  dplyr::select(Gene, average_phyloP)

library(dplyr)

consNTP_overall <- rnaseq %>%
  filter(Group == "Non-Triplets", 
         NMD == "Sensitive", 
         average_phyloP > 1,
         Gene != "NONE",  # Exclude "NONE"
         !grepl("\\.", Gene)  # Exclude Gene names with "."
  ) %>%
  dplyr::select(Gene, average_phyloP)







