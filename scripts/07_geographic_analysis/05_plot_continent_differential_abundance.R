# Purpose: Plot continent-enriched genera together with their prevalence and mean relative abundance.
# Panels: Figure 4A and Figure S8B.

library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggnewscale)

setwd("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/")
continent_list <- c("Asia", "North America", "Europe", "Africa")
mean_abun <- read.csv("differential_analysis_res/genus_p0.01_log10relative_abundance_mean_countinent.csv", row.names = 1) %>% rename(`North America` = North.America)
metadata <- read.csv("ID_meta_final.csv", sep = ",", header = T, row.names = 1) #%>% rename(Sample_ID = ID)
prevalence <- read.csv("differential_analysis_res/genus_p0.01_prevelance.csv", row.names = 1)
res <- read.csv("differential_analysis_res/plotdata_continent_diff_genus.csv", row.names = 1)

res_pre <- res %>%
  left_join(prevalence %>% rename_with(~c("Taxa", "Variable_pre", "Variable")), by = c("Taxa", "Variable")) %>%
  left_join(prevalence %>% rename_with(~c("Taxa", "Reference_pre", "Reference_variable")), by = c("Taxa", "Reference_variable")) %>%
  filter(Significance1 == "Sig") %>%
  select(Variable, Reference_variable, Taxa, Estimate, adjust.p.value, Variable_pre, Reference_pre, Significance2) %>%
  data.frame()  %>%
  filter(Significance2 == "Higher")

res2 <- res %>% filter(Significance2 == "Higher") %>%
  mutate(Label = paste0(Variable, " vs. ", Reference_variable)) %>%
  mutate(taxa = str_replace_all(Taxa, "g__", "")) %>%
  select(Label, taxa, p.value, adjust.p.value) %>%
  rbind(data.frame(Label = "North America vs. Europe", taxa = "Saccharomyces", p.value = NA, adjust.p.value = NA)) %>%
  mutate(significance = case_when(is.na(adjust.p.value) ~ 0, T ~ 1))
res3 <- res2 %>% select(Label, taxa, significance) %>%
  pivot_wider(names_from = "taxa", values_from = "significance") %>%
  pivot_longer(cols = unique(res2$taxa), names_to = "taxa", values_to = "significance") %>%
  mutate(significance = case_when(is.na(significance) | significance == 0 ~ "non-Sig", T ~ "Sig")) %>%
  mutate(temp = Label) %>%
  separate(col = temp, into = c("Continent", "Reference"), sep = " vs. ") %>%
  mutate(enrichment_label = case_when(significance == "Sig" ~ Continent, T ~ "non-Sig"))

pre1 <- prevalence %>% mutate(taxa = str_replace_all(Genus, "g__", "")) %>% filter(taxa %in% res2$taxa)
mean_abun_long <- mean_abun %>% rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = all_of(continent_list), names_to = "Continent", values_to = "Relative_abundance")
abun1 <- mean_abun_long %>% mutate(taxa = str_replace_all(Taxa, "g__", "")) %>% filter(taxa %in% res2$taxa)
abun2 <- abun1 %>% mutate(abundance = Relative_abundance) %>%
  select(taxa, abundance) %>%
  group_by(taxa) %>%
  summarise(abundance_mean = mean(abundance, na.rm = T), .groups = "drop") %>%
  ungroup %>%
  arrange(abundance_mean)

# order
taxa_list <- abun2$taxa
compare_list <- c(
  "Asia vs. North America", "Asia vs. Europe", "Asia vs. Africa",
  "North America vs. Asia", "North America vs. Europe", "North America vs. Africa",
  "Europe vs. Asia", "Europe vs. North America", "Europe vs. Africa",
  "Africa vs. Asia",  "Africa vs. North America", "Africa vs. Europe"
)
continent_list <- c("Asia",  "North America", "Europe", "Africa")
color_continent <- c("Asia" =          "#a41c36",
                     "North America" = "#08519d",
                     "Europe" =        "#682487",
                     "Africa" =        "#ef8b67")
color_list <- c("#a41c36", "#08519d", "#682487", "#ef8b67")


compare_summary <- table(res3$Label[res3$significance == "Sig"]) %>%
  data.frame() %>% rename(Label = Var1)
compare_summary$Continent <- c("Africa", "Africa","Africa",
                               "Asia", "Asia","Asia",
                               "Europe","Europe","Europe",
                               "North America","North America","North America")
compare_summary$Label <- factor(compare_summary$Label, levels = rev(compare_list))
p_compare_summary <- compare_summary %>%
  ggplot(aes(x = Freq, y = Label, fill = Continent)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(x = Freq, label = Freq),
            hjust = 0.5,
            size = 3,
            nudge_x = 1,
            color = "black") +
  xlim(0, 25)+
  labs(y = "", x = "Number of continent\nenriched genera") +
  scale_fill_manual(values = color_continent)+
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 1),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
    axis.title = element_text(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.length.y = unit(-1, "mm")
  )
p_compare_summary


# plot_compare
compare_order <- read.table("differential_analysis_res/continent_compare_order.txt", sep = "\t", header = T)
compare_order_long <- compare_order %>%
  pivot_longer(cols = c("Asia", "North.America", "Europe", "Africa"), names_to = "continent", values_to = "label") %>%
  filter(!is.na(label)) %>%
  mutate(label = as.factor(label), order = as.factor(order), continent = str_replace_all(continent, "North.America", "North America")) %>%
  arrange(order)
compare_order_long$order <- factor(compare_order_long$order, levels = rev(unique(compare_order_long$order)))
compare_order_long$continent <- factor(compare_order_long$continent, levels = continent_list)
compare_order_long <- compare_order_long %>% mutate(label2 = case_when(label == "1" ~ continent, T ~ "Reference continent"))

shade_df <- data.frame(
  xmin = c("Asia", "Europe"),
  xmax = c("Asia", "Europe"),
  ymin = -Inf,
  ymax = Inf
)

shade_df <- data.frame(continent = c("North America", "Africa"))
p_compare_order <-
  ggplot(data = compare_order_long, aes(x = continent, y = order)) +
  geom_point(aes(fill = label2), size = 4, shape = 21, alpha = 1)+
  scale_fill_manual(values = c(color_continent, "Reference continent" = "white"))+
  geom_tile(
    data = shade_df,
    aes(x = continent, y = 0.5, height = Inf, width = 1),
    fill = "grey", alpha = 0.2
  )+
  labs(title = "",
       x = "",
       y = "",
       fill = ""
  ) +
  theme_minimal()+
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.y = unit(0, "mm")
  )
p_compare_order

temp <- res %>%
  select(Variable, Taxa, Reference_variable, Estimate) %>%
  rename(Continent = Variable, Reference = Reference_variable, taxa = Taxa) %>%
  mutate(taxa = str_replace_all(taxa, "g__", ""))
res4 <- res3 %>%
  left_join(temp, by = c("Continent", "Reference", "taxa")) %>%
  mutate(estimate_for_plot = case_when(significance == "Sig" ~ Estimate, T ~ NA))
res4$Label <- factor(res4$Label, levels = rev(compare_list))
res4$taxa <- factor(res4$taxa, levels = taxa_list)


global_range <- range(res4$estimate_for_plot, na.rm = TRUE)
p_compare <- ggplot(res4, aes(x = taxa, y = Label)) +
  geom_tile(data = subset(res4, enrichment_label == "Africa"),
            aes(fill = estimate_for_plot), color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#f4e1c3", high = "#ef8b67", name = "Africa", limits = global_range, breaks = seq(0, 6, by = 2)) +

  new_scale_fill() +
  geom_tile(data = subset(res4, enrichment_label == "Europe"),
            aes(fill = estimate_for_plot), color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#ceabea", high = "#682487", name = "Europe", limits = global_range, breaks = seq(0, 6, by = 2)) +

  new_scale_fill() +
  geom_tile(data = subset(res4, enrichment_label == "North America"),
            aes(fill = estimate_for_plot), color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#a8cadf", high = "#08519d", name = "North America", limits = global_range, breaks = seq(0, 6, by = 2)) +

  new_scale_fill() +
  geom_tile(data = subset(res4, enrichment_label == "Asia"),
            aes(fill = estimate_for_plot), color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#f09ba0", high = "#a41c36", name = "Asia", limits = global_range, breaks = seq(0, 6, by = 2)) +

  new_scale_fill() +
  geom_tile(data = subset(res4, enrichment_label == "non-Sig"),
            fill = "#efefef",
            color = "white", linewidth = 0.5) +

  labs(title = "", x = "", y = "") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", face = "italic"),
    axis.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.y = unit(0, "mm")
  )
p_compare
p_compare_total <- p_compare_order + p_compare + p_compare_summary + plot_layout(nrow = 1, widths = c(1,  10, 2))
p_compare_total
ggsave(plot = p_compare_total, "plot/Differential_genus_compare2.0.pdf", width = 12, height = 5)







# plot prevalence & abundance
continent_list <- c("Asia", "North America", "Europe", "Africa")
pre1$Continent <- factor(pre1$Continent, levels = continent_list)
pre1$taxa <- factor(pre1$taxa, levels = taxa_list)

p_pre <- list()
for(i in 1:4){
  data0 <- pre1 %>% filter(Continent == continent_list[i]) %>%
    mutate(non_Pre = 1 - Pre) %>%
    pivot_longer(cols = c("Pre", "non_Pre"), names_to = "Type", values_to = "Prevalence")

  data0$taxa <- factor(data0$taxa, levels = taxa_list)
  p_pre[[i]] <- data0 %>%
    ggplot(aes(x = Prevalence, y = taxa, fill = Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Pre" = color_list[i],"non_Pre" = "#ebebeb"))+
    labs(title = "",
         x = continent_list[i],
         y = "",
         fill = ""
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 1))+
    theme_minimal()+
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
      axis.title = element_text(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(color = "black"),
      axis.ticks.length.y = unit(0, "mm")
    )
}
p_pre_merge <- p_pre[[1]] +  p_pre[[2]] + p_pre[[3]] + p_pre[[4]] + plot_layout(nrow = 1)
p_pre_merge

summary(abun1$Relative_abundance)
abun1$taxa <- factor(abun1$taxa, levels = taxa_list)
abun1$Continent <- factor(abun1$Continent, levels = continent_list)
p_abun <-
  ggplot(abun1, aes(x = Relative_abundance, y = taxa)) +
  geom_point(aes(fill = Continent), size = 2, shape = 21, alpha = 0.8)+
  scale_fill_manual(values = color_continent)+
  labs(title = "Relative abundance",
       x = "Mean(Log10(rel. abun.))",
       y = "",
       fill = ""
  ) +
  xlim(-6, 0)+
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line.x.bottom = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks.length.y = unit(0, "mm")
  )
p_abun

p_prev_abun <- p_abun + p_pre[[1]] + p_pre[[2]] + p_pre[[3]] + p_pre[[4]] + plot_layout(nrow = 1, widths = c(10, 2, 2, 2, 2))
p_prev_abun
ggsave(plot = p_prev_abun, "plot/Differential_genus_prevlance_abundance.pdf", width = 8, height = 8)
