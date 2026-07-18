# Purpose: Compare primary and sensitivity-analysis effect sizes and assemble the Figure S7 robustness panels.
# Panels: Figures S7A-C.

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)

# --- PERMANOVA sensitivity analysis for macrofungal filtering ---
setwd("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis")
perm <- read.table("sensitive_analysis/res_permanova_macrofungi_check.tsv", header = T)
results <- perm %>%
  group_by(diversity_index) %>%
  mutate(
    p_anova_adj = p.adjust(p.value, method = "fdr")
  ) %>%
  ungroup()

effect_size_meta <- read.csv("sensitive_analysis/effect_size_meta.csv")
permanova_effect_meta_macrofungi <- merge(results, effect_size_meta, by = "factor", all.x = T) %>%
  filter(to_analysis == 1) %>%
  mutate(sig_label = case_when(p_anova_adj > 0.01 ~ "ns", T ~ NA))

permanova_effect_meta_agar <- permanova_effect_meta_macrofungi %>%
  filter(diversity_index == "agar") %>%
  arrange(R2) %>%
  mutate(name = factor(name, levels = unique(name)))

permanova_effect_meta_list <- permanova_effect_meta_macrofungi %>%
  filter(diversity_index == "list") %>%
  arrange(R2) %>%
  mutate(name = factor(name, levels = unique(name)))

type_colors <- c("Demography" = "#001427",
                 "Geography" = "#708d81",
                 "Sequencing" = "#f4d58d",
                 "DNA_extract" = "#274c77",
                 "PCR" = "#8d0801",
                 "Sample_type" = "#999999")

res_permanova <- read.csv("sensitive_analysis/res_PERMANOVA.csv") # all data result
merged_data_agar <- permanova_effect_meta_agar %>%
  select(factor, R2, name, type) %>%
  unique() %>%
  left_join(res_permanova %>% rename(total_R2 = R2), by = "factor")
pearson_cor_agar <- cor.test(merged_data_agar$total_R2, merged_data_agar$R2, method = "pearson")
p_value_agar <- pearson_cor_agar$p.value; cor_coef_agar <- pearson_cor_agar$estimate
p_cor_agar <- ggplot(merged_data_agar, aes(y = total_R2*100, x = R2*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (R"^2, ", all data)")),
       x = expression(paste("Effect size (R"^2, ", micro_agar)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 6,
           label = paste0(
             "atop(\"Pearson's r\" == ", format(cor_coef_agar, digits=6),
             ", italic(P) == ", format(p_value_agar, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.position = "none"
  )
p_cor_agar


merged_data_list <- permanova_effect_meta_list %>%
  select(factor, R2, name, type) %>%
  unique() %>%
  left_join(res_permanova %>% rename(total_R2 = R2), by = "factor")
pearson_cor_list <- cor.test(merged_data_list$total_R2, merged_data_list$R2, method = "pearson")
p_value_list <- pearson_cor_list$p.value; cor_coef_list <- pearson_cor_list$estimate
p_cor_list <- ggplot(merged_data_list, aes(y = total_R2*100, x = R2*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (R"^2, ", all data)")),
       x = expression(paste("Effect size (R"^2, ", micro_list)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 6,
           label = paste0(
             "atop(\"Pearson's r\" == ", format(cor_coef_list, digits=6),
             ", italic(P) == ", format(p_value_list, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.position = "none"
  )
p_cor_list



# --- Alpha-diversity effect-size sensitivity analysis for macrofungal filtering ---
final_result <- read.csv("sensitive_analysis/effect_size_diversity_lm_macro.csv")
res_for_cor <- final_result %>%
  select(diversity_index, eta_squared, factor) %>%
  pivot_wider(names_from = diversity_index, values_from = eta_squared) %>%
  left_join(effect_size_meta, by = "factor") %>%
  filter(to_analysis == 1)

pearson_cor_lm_list <- cor.test(res_for_cor$Shannon_micro_dual, res_for_cor$Shannon_micro_list, method = "pearson")
p_value_lm_list <- pearson_cor_lm_list$p.value; cor_coef_lm_list <- pearson_cor_lm_list$estimate
p_cor_lm_list <- ggplot(res_for_cor, aes(y = Shannon_micro_dual*100, x = Shannon_micro_list*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (" * eta^2 * ", all data)")),
       x = expression(paste("Effect size (" * eta^2 * ", micro_list)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 15,
           label = paste0(
             "atop(\"Pearson's r\" == ", format(cor_coef_lm_list, digits=6),
             ", italic(P) == ", format(p_value_lm_list, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.position = "none"
  )
p_cor_lm_list


pearson_cor_lm_agar <- cor.test(res_for_cor$Shannon_micro_dual, res_for_cor$Shannon_micro_agar, method = "pearson")
p_value_lm_agar <- pearson_cor_lm_agar$p.value; cor_coef_lm_agar <- pearson_cor_lm_agar$estimate
p_cor_lm_agar <- ggplot(res_for_cor, aes(y = Shannon_micro_dual*100, x = Shannon_micro_agar*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (" * eta^2 * ", all data)")),
       x = expression(paste("Effect size (" * eta^2 * ", micro_agar)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 15,
           label = paste0(
             "atop(\"Pearson's r\" == ", format(cor_coef_lm_agar, digits=6),
             ", italic(P) == ", format(p_value_lm_agar, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.position = "none"
  )
p_cor_lm_agar









# --- Alpha-diversity effect sizes in balanced samples ---
results <- read.csv("sensitive_analysis/effect_size_diversity_lm_subsample.csv")
effect_size_meta <- read.csv("sensitive_analysis/effect_size_meta.csv")
lm_effect_meta <- merge(results, effect_size_meta, by = "factor", all.x = T) %>%
  filter(to_analysis == 1) %>%
  mutate(sig_label = case_when(p_anova_adj > 0.01 ~ "ns", T ~ NA)) %>%
  group_by(name) %>%
  mutate(eta_squared_mean = mean(eta_squared)) %>%
  ungroup() %>%
  arrange(eta_squared_mean) %>%
  mutate(name = factor(name, levels = unique(name)))

res_lm <- read.csv("sensitive_analysis/effect_size_diversity_lm.csv") # all data
merged_data_lm <- lm_effect_meta %>%
  select(factor, eta_squared, name, type) %>%
  unique() %>%
  left_join(res_lm %>% filter(div_index == "Shannon") %>% select(factor, eta_squared) %>% rename(total_eta_squared = eta_squared), by = "factor")

pearson_cor <- cor.test(merged_data_lm$total_eta_squared, merged_data_lm$eta_squared, method = "pearson")
p_value <- pearson_cor$p.value; cor_coef <- pearson_cor$estimate

merged_data_lm <- merged_data_lm %>% group_by(name) %>% arrange(eta_squared) %>% mutate(name_label = if_else(row_number() == 1, name, NA)) %>% ungroup
p_cor_lm_subsample <- ggplot(merged_data_lm, aes(y = total_eta_squared*100, x = eta_squared*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (" * eta^2 * ", all data)")),
       x = expression(paste("Effect size (" * eta^2 * ", balanced samples)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 15,
           label = paste0(
             "atop(\"Pearson's r\" == ", format(cor_coef, digits=6),
             ", italic(P) == ", format(p_value, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.key.size = unit(0.15, "in")
  )
p_cor_lm_subsample


# --- PERMANOVA effect sizes in balanced samples ---
perm <- read.table("sensitive_analysis/res_permanova_subsample.tsv", header = T)
results <- perm %>%
  arrange(factor, subsample) %>%
  group_by(subsample) %>%
  mutate(
    p_anova_adj = p.adjust(p.value, method = "fdr")
  ) %>%
  ungroup()

effect_size_meta <- read.csv("sensitive_analysis/effect_size_meta.csv")
permanova_effect_meta <- merge(results, effect_size_meta, by = "factor", all.x = T) %>%
  filter(to_analysis == 1) %>%
  mutate(sig_label = case_when(p_anova_adj > 0.01 ~ "ns", T ~ NA)) %>%
  group_by(name) %>%
  mutate(R2_mean = mean(R2)) %>%
  ungroup() %>%
  arrange(R2_mean) %>%
  mutate(name = factor(name, levels = unique(name)))

res_permanova <- read.csv("sensitive_analysis/res_PERMANOVA.csv") # all samples result
merged_data_perm <- permanova_effect_meta %>%
  select(factor, R2, name, type) %>%
  unique() %>%
  left_join(res_permanova %>% rename(total_R2 = R2), by = "factor")

pearson_cor_perm <- cor.test(merged_data_perm$total_R2, merged_data_perm$R2, method = "pearson")
p_value_perm <- pearson_cor_perm$p.value; cor_coef_perm <- pearson_cor_perm$estimate


merged_data_perm <- merged_data_perm %>% group_by(name) %>% arrange(R2) %>% mutate(name_label = if_else(row_number() == 1, name, NA)) %>% ungroup
p_cor_perm_subsample <- ggplot(merged_data_perm, aes(y = total_R2*100, x = R2*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (R"^2, ", all data)")),
       x = expression(paste("Effect size (R"^2, ", balanced samples)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 6,
           label = paste0(
             "atop(\"Pearson's r\" == ", format(cor_coef_perm, digits=6),
             ", italic(P) == ", format(p_value_perm, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.position = "none"
  )
p_cor_perm_subsample



# --- Geographic PERMANOVA adjusted for individual technical factors ---
res_adjusttech <- read.table("sensitive_analysis/res_permsnova_adjusttech.tsv", header = T) %>%
  group_by(factor, adj_factor) %>%
  slice(2) %>%
  ungroup %>%
  rename(p.value = Pr..F.) %>%
  select(factor, adj_factor, R2, p.value) %>%
  arrange(factor) %>% data.frame()
res_adjusttech

res_allsample <- read.csv("sensitive_analysis/res_PERMANOVA.csv", header = T) %>% rename(p.value = p_value)
factor_list <- c("Continent", "Country", "Region")
merge_data_allsample <- res_allsample %>%
  mutate(label = "Unadjusted model") %>%
  bind_rows(res_adjusttech %>% mutate(label = "Technical factor adjusted model")) %>%
  filter(factor %in% factor_list)
merge_data_allsample
plot_data <- merge_data_allsample %>%
  group_by(factor, label) %>%
  summarise(
    Mean_R2 = mean(R2),
    SD_R2 = sd(R2),
    .groups = 'drop'
  ) %>%
  mutate(
    y = Mean_R2 * 100,
    lower_ci = y - (SD_R2 * 100),
    upper_ci = y + (SD_R2 * 100)
  )
plot_data

plot_data$label <- factor(plot_data$label, levels = c("Unadjusted model", "Technical factor adjusted model"))
plot_data$factor <- factor(plot_data$factor, levels = c("Region", "Country", "Continent"))
merge_data_allsample$label <- factor(merge_data_allsample$label, levels = c("Unadjusted model", "Technical factor adjusted model"))
p_allsample_techadj <- ggplot() +
  geom_bar(data = plot_data,
           aes(x = factor, y = y, fill = label, group = label),
           stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.8, alpha = 0.9) +
  geom_jitter(data = merge_data_allsample,
              aes(x = factor, y = R2 * 100,  group = label),
              color = "black", fill = "#708d81", shape  = 21,
              position = position_jitterdodge(jitter.width = 0.2,
                                              jitter.height = 0,
                                              dodge.width = 0.8),
              size = 1.5,
              alpha = 0.25) +
  geom_text(data = plot_data,
            aes(x = factor, y = y, label = paste0(round(y, 1), "%"), group = label),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 10),  breaks = seq(0, 10, by = 2)) +
  scale_fill_manual(values = c("Unadjusted model" = "#708d81",
                               "Technical factor adjusted model" = "#cfdfc5")) +
  scale_color_manual(values = c("Unadjusted model" = "#708d81",
                                "Technical factor adjusted model" = "#cfdfc5")) +
  labs(y = expression(paste("Effect size (R"^2, ", all data)")), x = "", fill = NULL, color = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.5, 0.9),
        legend.key.size = unit(0.15, "in"))
p_allsample_techadj


# --- PERMANOVA sensitivity analysis after removing Unknown ---
perm <- read.table("sensitive_analysis/res_permanova_dropunknown.tsv", header = T)
results <- perm %>%
  group_by(diversity_index) %>%
  mutate(
    p_anova_adj = p.adjust(p.value, method = "fdr")
  ) %>%
  ungroup()

effect_size_meta <- read.csv("sensitive_analysis/effect_size_meta.csv")
permanova_effect_meta_dropunknown <- merge(results, effect_size_meta, by = "factor", all.x = T) %>%
  filter(to_analysis == 1) %>%
  mutate(sig_label = case_when(p_anova_adj > 0.01 ~ "ns", T ~ NA)) %>%
  group_by(name) %>%
  mutate(R2_mean = mean(R2)) %>%
  ungroup() %>%
  arrange(R2_mean) %>%
  mutate(name = factor(name, levels = unique(name)))

res_permanova <- read.csv("sensitive_analysis/res_PERMANOVA.csv") # all data result
merged_data_dropunknown <- permanova_effect_meta_dropunknown %>%
  select(factor, R2, name, type) %>%
  unique() %>%
  left_join(res_permanova %>% rename(total_R2 = R2), by = "factor")
head(merged_data_dropunknown)

pearson_cor <- cor.test(merged_data_dropunknown$total_R2, merged_data_dropunknown$R2, method = "pearson")
p_value <- pearson_cor$p.value; cor_coef <- pearson_cor$estimate
p_cor_dropunknown <- ggplot(merged_data_dropunknown, aes(y = total_R2*100, x = R2*100, fill = type)) +
  geom_point(alpha = 0.9, size = 4, shape = 21, color = "black")+
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  geom_smooth(method = "lm", se = T, color = "#e94f37", linewidth = 0.75, formula = y ~ x, fullrange = TRUE, fill = "grey", alpha = 0.5) +
  scale_fill_manual(values = type_colors) +
  labs(y = expression(paste("Effect size (R"^2, ", all data)")),
       x = expression(paste("Effect size (R"^2, ", dropunknown)")),
       fill = NULL
  ) +
  annotate("text",
           x = 0, y = 6,
           label = paste0(
             "atop(Pearson_rho == ", format(cor_coef, digits=6),
             ", italic(P) == ", format(p_value, digits=2), ")"
           ),
           hjust = 0, vjust = 1, size = 3,
           parse = TRUE) +
  theme_minimal()+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.x = unit(1, "mm"),
    axis.ticks.length.y = unit(1, "mm"),
    legend.key.size = unit(0.15, "in")
  )
p_cor_dropunknown


# --- Assemble the complete sensitivity figure ---
p_lm <- p_cor_lm_agar + p_cor_lm_list + p_cor_lm_subsample + p_allsample_techadj + plot_layout(nrow = 1)
p_permanova <- p_cor_agar + p_cor_list + p_cor_perm_subsample + p_cor_dropunknown + plot_layout(nrow = 1)
plot_all <- p_lm / p_permanova
plot_all
ggsave(plot = plot_all, "sensitive_analysis.pdf", width = 12, height = 7)
