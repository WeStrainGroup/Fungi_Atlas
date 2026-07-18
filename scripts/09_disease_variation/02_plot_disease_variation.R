# Purpose: Plot cohort diversity, PCoA, differential abundance, prediction performance, and the composite disease figure.
# Panels: Figure 6A-C; Figure S9A-C.

library(tidyverse)
library(ggplot2)
library(patchwork)
library(hugmycoa)
library(phyloseq)


### plot sample number
setwd("/data/huangkailang/project/fungi_ITS/case_control_metadata_clean")

sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)
cohort_order
sample_number$cohort <- factor(sample_number$cohort, levels = cohort_order)
y_levels <- unique(sample_number$cohort)
rect_data <- data.frame(y_idx = seq_along(y_levels)) %>% filter(y_idx %% 2 == 0)

p_sample_number <-
  sample_number %>%
  ggplot(aes(y = cohort, x = n_sample, fill = group)) +
  geom_rect(data = rect_data,
            aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5, xmin = -Inf, xmax = Inf),
            fill = "gray95",
            inherit.aes = FALSE) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Control" = "#bebdbd", "Case" = "#96bccb")) +
    labs(
    x = "Number of samples",
    y = "",
    fill = ""
  ) +
  theme_minimal() +
  theme(
        legend.position = c(0.8, 0.1),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        axis.title.x = element_text(size = 8)
      )
p_sample_number



### plot alpha diversity wilcox test
result_alpha_all <- read.csv("table/res_alpha_diversity_wilcox_test.csv")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)

plot <- list()
for(index in c("Shannon", "Simpson", "Observed")) {
    df0 <- result_alpha_all %>% filter(alpha_index == index) %>%
        mutate(cohort = factor(cohort, levels = cohort_order))  %>%
        arrange(cohort)
    y_levels <- unique(df0$cohort)
    rect_data <- data.frame(y_idx = seq_along(y_levels)) %>%
                 filter(y_idx %% 2 == 0)

    plot[[index]] <- df0 %>%
      ggplot(aes(y = cohort, x = Estimate, fill = alpha_index)) +
      geom_rect(data = rect_data,
                aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                    xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
      geom_col(color = "white", linewidth = 0.3) +
      scale_fill_manual(values = c("Shannon" = "#d88a7a", "Simpson" = "#d88a7a", "Observed" = "#d88a7a")) +
      geom_text(
        aes(
          label = significance,
          hjust = ifelse(Estimate >= 0, -0.2, 1.2)
          ),
        size = 4, vjust = 0.8
        ) +
      labs(x = if(index == "Observed") paste0("Estimated difference \n(", index, " feature)") else paste0("Estimated difference \n(", index, " index)"),
           y = "") +
      scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 8)
      )
}
p_alpha <- plot[["Observed"]] + plot[["Shannon"]] + plot[["Simpson"]] +  plot_layout(guides = "collect", nrow = 1)



### plot beta diversity PERMANOVA test (genus level)
result_adonis_robust_genus <- read.csv("table/res_beta_diversity_adonis_test_robust_genus.csv")
y_levels <- unique(result_adonis_robust_genus$cohort)
rect_data <- data.frame(y_idx = seq_along(y_levels)) %>% filter(y_idx %% 2 == 0)
p_adonis_robust_genus <- result_adonis_robust_genus %>%
      ggplot(aes(y = cohort, x = adonis_R2 * 100)) +
      geom_rect(data = rect_data,
                aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                    xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
      geom_col(fill = "#95c1b9", color = "white", linewidth = 0.3) +
      geom_text(
        aes(
          label = significance,
          hjust = -0.2
        ),
        size = 4,
        vjust = 0.8
      ) +
      labs(title = "", x = paste0("Adonis R², %\n(Genus; Robust Aitchison distance)"), y = "") +
      xlim(0, 70)+

      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = if(index == "Shannon") element_text(color = "black", size = 8) else element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = if(index == "Shannon") element_line(color = "black") else element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 8)
      )
p_adonis_robust_genus



### plot beta diversity PERMANOVA test (SH level)
result_adonis_robust_sh <- read.csv("table/res_beta_diversity_adonis_test_robust_sh.csv")
y_levels <- unique(result_adonis_robust_sh$cohort)
rect_data <- data.frame(y_idx = seq_along(y_levels)) %>% filter(y_idx %% 2 == 0)
p_adonis_robust_sh <- result_adonis_robust_sh %>%
      ggplot(aes(y = cohort, x = adonis_R2 * 100)) +
      geom_rect(data = rect_data,
                aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                    xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
      geom_col(fill = "#95c1b9", color = "white", linewidth = 0.3) +
      geom_text(
        aes(
          label = significance,
          hjust = -0.2
        ),
        size = 4,
        vjust = 0.8
      ) +
      labs(title = "", x = paste0("Adonis R², %\n(SH; Robust Aitchison distance)"), y = "") +
      xlim(0, 70)+

      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = if(index == "Shannon") element_text(color = "black", size = 8) else element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = if(index == "Shannon") element_line(color = "black") else element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 8)
      )
p_adonis_robust_sh



### plot beta diversity pcoa （genus)
res_beta_diversity_pcoa_genus <- read.csv("table/res_beta_diversity_pcoa_genus.csv")
cohort_list <- cohort_order
pcoa_genus <- list()
for(cohort0 in cohort_list) {
  df0 <- res_beta_diversity_pcoa_genus %>% filter(cohort == cohort0)
adonis2_label <- paste0("Adonis R² = ", round(df0$adonis_R2[1] * 100, 2), " %  ",
                        "P ", ifelse(df0$adonis_p_value[1] < 0.001, "< 0.001", paste0("= ", round(df0$adonis_p_value[1], 3))))
pcoa_genus[[cohort0]] <- df0 %>%
  ggplot(aes(x = PCoA1, y = PCoA2, fill = group, color = group)) +
  geom_point(size = 2.5, shape = 16, alpha = 0.9) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.1,
    level = 0.95,
    type = "t",
    linewidth = 0.5) +
  scale_color_manual(values = c("Control" = "#6a91b8", "Case" = "#b3656e")) +
  scale_fill_manual(values = c("Control" = "#6a91b8", "Case" = "#b3656e")) +
  labs(x = paste0("PCoA1 (", df0$pcoa_eig1[1], "%)"),
       y = paste0("PCoA2 (", df0$pcoa_eig2[1], "%)"),
       color="", fill="", title = paste0(cohort0, "\n", adonis2_label))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.9, 0.1),
        legend.key.size = unit(0.15, "in")) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))
}
p_pcoa_genus <- wrap_plots(pcoa_genus[1:23], ncol = 5, guides = "collect") &
     theme(legend.position = "right")
p_pcoa_genus
ggsave(plot = p_pcoa_genus, "plot/pcoa_plot_genus.pdf", width = 12, height = 12)



### plot beta diversity pcoa （sh)
res_beta_diversity_pcoa_sh <- read.csv("table/res_beta_diversity_pcoa_sh.csv")
cohort_list <- cohort_order
pcoa_sh <- list()
for(cohort0 in cohort_list) {
  df0 <- res_beta_diversity_pcoa_sh %>% filter(cohort == cohort0)
adonis2_label <- paste0("Adonis R² = ", round(df0$adonis_R2[1] * 100, 2), " %  ",
                        "P ", ifelse(df0$adonis_p_value[1] < 0.001, "< 0.001", paste0("= ", round(df0$adonis_p_value[1], 3))))
pcoa_sh[[cohort0]] <- df0 %>%
  ggplot(aes(x = PCoA1, y = PCoA2, fill = group, color = group)) +
  geom_point(size = 2.5, shape = 16, alpha = 0.9) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.1,
    level = 0.95,
    type = "t",
    linewidth = 0.5) +
  scale_color_manual(values = c("Control" = "#6a91b8", "Case" = "#b3656e")) +
  scale_fill_manual(values = c("Control" = "#6a91b8", "Case" = "#b3656e")) +
  labs(x = paste0("PCoA1 (", df0$pcoa_eig1[1], "%)"),
       y = paste0("PCoA2 (", df0$pcoa_eig2[1], "%)"),
       color="", fill="", title = paste0(cohort0, "\n", adonis2_label))+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.9, 0.1),
        legend.key.size = unit(0.15, "in")) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))
}
p_pcoa_sh <- wrap_plots(pcoa_sh[1:23], ncol = 5, guides = "collect") &
     theme(legend.position = "right")
p_pcoa_sh
ggsave(plot = p_pcoa_sh, "plot/pcoa_plot_sh.pdf", width = 12, height = 12)



### plot differential analysis (genus)
library(ggnewscale)
setwd("/data/huangkailang/project/fungi_ITS/case_control_metadata_clean")
res_diff_genus <- read.csv("table/ancombc2_results_genus.csv")

res_sum <- res_diff_genus %>%
  filter(abs(log2FC) > 1 & q_value < 0.05)  %>%
  mutate(diff_type = ifelse(log2FC > 0, "Case_enriched", "Control_enriched")) %>%
  group_by(feature_id, diff_type)  %>%
  summarise(count = n(), .groups = "drop") %>%
  ungroup() %>%
  group_by(feature_id) %>%
  mutate(total_count = sum(count),
         diff_type_count = n()) %>%
  mutate(
    case_count = sum(count[diff_type == "Case_enriched"], na.rm = TRUE),
    control_count = sum(count[diff_type == "Control_enriched"], na.rm = TRUE),
    feature_type = case_when(
      diff_type_count == 1 ~ as.character(diff_type),
      case_count > control_count ~ "Case_enriched",
      case_count < control_count ~ "Control_enriched",
      TRUE ~ "Cmixed_enriched"
    )
  ) %>%
  ungroup() %>%
  mutate(total_count1 = ifelse(feature_type == "Case_enriched", total_count, -total_count),
         menus_count = case_count - control_count) %>%
  mutate(diff_type = factor(diff_type, levels = c("Case_enriched", "Control_enriched")))  %>%
  data.frame() %>%
  arrange(-total_count1) %>%
  mutate(feature_label = str_remove_all(feature_id, "g__"))


feature_order_table <- res_sum %>%
  select(feature_label, feature_type, total_count1, menus_count) %>%
  unique() %>%
  mutate(total_count1 = as.numeric(total_count1), menus_count = as.numeric(menus_count)) %>%
  arrange(feature_type, desc(total_count1), desc(menus_count))
head(feature_order_table)

res_sum$feature_label <- factor(res_sum$feature_label, levels = feature_order_table$feature_label)
cohort_sum <- res_diff_genus %>%
  group_by(cohort) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  select(cohort,count) %>%
  unique() %>%
  arrange(-count) %>%
  mutate(cohort = factor(cohort, levels = unique(cohort)))

feature_levels <- levels(res_sum$feature_label)
n_feat <- length(feature_levels)
x_pos <- seq_len(n_feat)
shade_positions <- x_pos[x_pos %% 2 == 1]
shade_df <- data.frame(
  xmin = shade_positions - 0.5,
  xmax = shade_positions + 0.5,
  ymin = -Inf,
  ymax = Inf
)

p_sig_feature_summary <- ggplot(res_sum, aes(x = feature_label, y = count, fill = diff_type)) +
  geom_rect(data = shade_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#f6f6f6", alpha = 1, inherit.aes = FALSE) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3, alpha = 1) +
  labs(x = "", y = "", fill = "Difference Type") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", size = 6, face = "italic"),
    axis.text.y = element_text(color = "black", size = 8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black")
  ) +
  scale_fill_manual(values = c("Case_enriched" = "#b3656e", "Control_enriched" = "#6a91b8"))
p_sig_feature_summary


res_diff_genus <- read.csv("table/ancombc2_results_genus.csv")
sig_features <- res_diff_genus %>%
  filter(abs(log2FC) > 1 & q_value < 0.05) %>%
  pull(feature_id) %>%
  unique()
res_diff_genus1 <- res_diff_genus %>%
  filter(feature_id %in% sig_features) %>%
  mutate(significance = ifelse(!is.na(log2FC) & abs(log2FC) > 1 & q_value < 0.05, "significant", "non-significant")) %>%
  mutate(cohort = factor(cohort, levels = unique(cohort))) %>%
  mutate(significance2 = case_when(
    significance == "significant" & log2FC > 0 & q_value < 0.05 ~ "Case_enriched",
    significance == "significant" & log2FC < 0 & q_value < 0.05 ~ "Control_enriched",
    TRUE ~ "Non-significant"
  ))

add_data <- data.frame()
for(cohort0 in unique(res_diff_genus1$cohort)) {
  feature_id0 <- res_diff_genus1 %>% filter(cohort == cohort0) %>% pull(feature_id) %>% unique()
  add <- setdiff(unique(res_diff_genus1$feature_id), feature_id0)
  add_data0 <- data.frame(cohort = cohort0, feature_id = add)
  add_data <- rbind(add_data, add_data0)
}
add_data$significance2 <- "Not_detected"
res_plot <- res_diff_genus1 %>%
  bind_rows(add_data) %>%
  mutate(feature_label = str_remove_all(feature_id, "g__")) %>%
  mutate(feature_label = factor(feature_label, levels = feature_order_table$feature_label)) %>%
  mutate(cohort = factor(cohort, levels = cohort_order))

p_sig_feature <- ggplot(res_plot, aes(y = cohort, x = feature_label)) +
  geom_tile(data = subset(res_plot, significance2 == "Non-significant"), fill = "#dddddd", color = "white", linewidth = 0.2) +
  new_scale_fill() +
  geom_tile(data = subset(res_plot, significance2 == "Not_detected"), fill = "#f6f6f6", color = "white", linewidth = 0.2) +
  new_scale_fill() +
  geom_tile(data = subset(res_plot, significance2 %in% c("Case_enriched", "Control_enriched")), aes(fill = log2FC), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#336190", mid = "#f6f6f6", high = "#992734", midpoint = 0,
                       limits = c(-max(abs(res_plot$log2FC), na.rm = TRUE), max(abs(res_plot$log2FC), na.rm = TRUE))) +
  scale_x_discrete(limits = levels(res_plot$feature_label)) +
  labs(x = "", y = "", fill = "Log2(FoldChange)",  title = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(color = "black", size = 8),
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_line(color = "black"),
    legend.text = element_text(size = 8)
  )
p_sig_feature
plot2 <- p_sig_feature / p_sig_feature_summary + plot_layout(heights = c(4, 1), guides = "collect")



### plot differential analysis (sh)
library(ggnewscale)
res_diff_sh <- read.csv("table/ancombc2_results_SH.csv")
sh_label_table <- read.csv("table/sh_label_table.csv")
res_sum <- res_diff_sh %>%
  filter(abs(log2FC) > 1 & q_value < 0.05)  %>%
  mutate(diff_type = ifelse(log2FC > 0, "Case_enriched", "Control_enriched")) %>%
  group_by(feature_id, diff_type)  %>%
  summarise(count = n(), .groups = "drop") %>%
  ungroup() %>%
  group_by(feature_id) %>%
  mutate(total_count = sum(count),
         diff_type_count = n()) %>%
  mutate(
    case_count = sum(count[diff_type == "Case_enriched"], na.rm = TRUE),
    control_count = sum(count[diff_type == "Control_enriched"], na.rm = TRUE),
    feature_type = case_when(
      diff_type_count == 1 ~ as.character(diff_type),
      case_count > control_count ~ "Case_enriched",
      case_count < control_count ~ "Control_enriched",
      TRUE ~ "Cmixed_enriched"
    )
  ) %>%
  ungroup() %>%
  mutate(total_count1 = ifelse(feature_type == "Case_enriched", total_count, -total_count),
         menus_count = case_count - control_count) %>%
  mutate(diff_type = factor(diff_type, levels = c("Case_enriched", "Control_enriched")))  %>%
  data.frame() %>%
  arrange(-total_count1) %>%
  mutate(feature_label = sh_label_table$sh_label[match(feature_id, sh_label_table$SH)])
head(res_sum)

feature_order_table <- res_sum %>%
  select(feature_label, feature_type, total_count1, menus_count) %>%
  unique() %>%
  mutate(total_count1 = as.numeric(total_count1), menus_count = as.numeric(menus_count)) %>%
  arrange(feature_type, desc(total_count1), desc(menus_count))
head(feature_order_table)

res_sum$feature_label <- factor(res_sum$feature_label, levels = feature_order_table$feature_label)
cohort_sum <- res_diff_sh %>%
  group_by(cohort) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  select(cohort,count) %>%
  unique() %>%
  arrange(-count) %>%
  mutate(cohort = factor(cohort, levels = unique(cohort)))


feature_levels <- levels(res_sum$feature_label)
n_feat <- length(feature_levels)
y_pos <- seq_len(n_feat)
shade_positions <- y_pos[y_pos %% 2 == 1]
shade_df <- data.frame(
  ymin = shade_positions - 0.5,
  ymax = shade_positions + 0.5,
  xmin = -Inf,
  xmax = Inf
)

res_sum <- res_sum %>%
  mutate(
    species = sub("^(.*)\\s+(ITS\\d+_SH\\d+)$", "\\1", feature_label),
    seq_id = sub("^(.*)\\s+(ITS\\d+_SH\\d+)$", "\\2", feature_label),
    formatted_label = paste0("italic('", species, "')~", seq_id)
  )

p_sig_feature_summary <- ggplot(res_sum, aes(y = feature_label, x = count, fill = diff_type)) +
  geom_rect(data = shade_df,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "#f6f6f6", alpha = 1, inherit.aes = FALSE) +
  geom_col(position = "stack", color = "white", linewidth = 0.3, alpha = 1) +
  labs(x = "", y = "", fill = "Difference Type") +
  scale_y_discrete(labels = function(breaks) {
    idx <- match(breaks, res_sum$feature_label)
    parse(text = res_sum$formatted_label[idx])
  }) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = "black", size = 5, face = "plain"),
    axis.text.x = element_text(color = "black", size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black")
  ) +
  scale_fill_manual(values = c("Case_enriched" = "#b3656e", "Control_enriched" = "#6a91b8"))
p_sig_feature_summary

res_diff_sh <- read.csv("table/ancombc2_results_SH.csv")
sig_features <- res_diff_sh %>%
  filter(abs(log2FC) > 1 & q_value < 0.05) %>%
  pull(feature_id) %>%
  unique()
res_diff_sh1 <- res_diff_sh %>%
  filter(feature_id %in% sig_features) %>%
  mutate(significance = ifelse(!is.na(log2FC) & abs(log2FC) > 1 & q_value < 0.05, "significant", "non-significant")) %>%
  mutate(cohort = factor(cohort, levels = unique(cohort))) %>%
  mutate(significance2 = case_when(
    significance == "significant" & log2FC > 0 & q_value < 0.05 ~ "Case_enriched",
    significance == "significant" & log2FC < 0 & q_value < 0.05 ~ "Control_enriched",
    TRUE ~ "Non-significant"
  ))

add_data <- data.frame()
for(cohort0 in unique(res_diff_sh1$cohort)) {
  feature_id0 <- res_diff_sh1 %>% filter(cohort == cohort0) %>% pull(feature_id) %>% unique()
  add <- setdiff(unique(res_diff_sh1$feature_id), feature_id0)
  add_data0 <- data.frame(cohort = cohort0, feature_id = add)
  add_data <- rbind(add_data, add_data0)
}
add_data$significance2 <- "Not_detected"
res_plot <- res_diff_sh1 %>% bind_rows(add_data) %>%
  mutate(feature_label = sh_label_table$sh_label[match(feature_id, sh_label_table$SH)])

res_plot$feature_label <- factor(res_plot$feature_label, levels = feature_order_table$feature_label)
res_plot$cohort <- factor(res_plot$cohort, levels = cohort_order)

table(res_plot$significance2)

p_sig_feature <- ggplot(res_plot, aes(x = cohort, y = feature_label)) +
  geom_tile(data = subset(res_plot, significance2 == "Non-significant"), fill = "#dddddd", color = "white", linewidth = 0.2) +
  new_scale_fill() +
  geom_tile(data = subset(res_plot, significance2 == "Not_detected"), fill = "#f6f6f6", color = "white", linewidth = 0.2) +
  new_scale_fill() +
  geom_tile(data = subset(res_plot, significance2 %in% c("Case_enriched", "Control_enriched")), aes(fill = log2FC), color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#336190", mid = "#f6f6f6", high = "#992734", midpoint = 0,
                       limits = c(-max(abs(res_plot$log2FC), na.rm = TRUE), max(abs(res_plot$log2FC), na.rm = TRUE))) +
  scale_y_discrete(limits = levels(res_plot$feature_label)) +
  labs(x = "", y = "", fill = "Log2(FoldChange)",  title = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", size = 6),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_line(color = "black"),
    legend.text = element_text(size = 6)
  )
p_sig_feature

plot_diff_sh <- (p_sig_feature_summary | p_sig_feature) + plot_layout(widths = c(1, 4), guides = "collect")
plot_diff_sh



### plot genus level prediction
res_pre_genus <- read.csv("table/res_prediction_genus.csv")
res_pre_genus <- res_pre_genus %>%
  filter(model %in% c("lasso", "rf")) %>%
  mutate(mechine_model  = str_replace_all(model, c( "svm" = "SVM", "rf" = "Random forest", "lasso" = "Lasso")))
res_pre_genus$cohort <- factor(res_pre_genus$cohort, levels = cohort_order)
res_pre_genus$mechine_model = factor(res_pre_genus$mechine_model, levels = rev(c("Random forest", "SVM", "Lasso")))

y_levels <- unique(res_pre_genus$cohort)
rect_data <- data.frame(y_idx = seq_along(y_levels)) %>% filter(y_idx %% 2 == 0)
p_n_feature_genus <- res_pre_genus %>% select(cohort, n_features) %>% unique() %>%
  ggplot(aes(y = cohort, x = n_features)) +
  geom_rect(data = rect_data,
            aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
  geom_col(fill = "#b2b6c1", color = "white", linewidth = 0.3) +
  labs(title = "",
       y = "",
       x = "Number of feature\n(Genus level)"
       ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "black")
        )
p_n_feature_genus

p_performance_genus <- res_pre_genus %>%
  ggplot(aes(y = cohort, x = auc, fill = mechine_model)) +
  geom_rect(data = rect_data,
            aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(xmin = low, xmax = high),
                position = position_dodge(width = 0.8),
                width = 0, linewidth = 0.1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#212121", lwd = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("SVM" = "#8cb8dc",
                               "Random forest" = "#bfdcdc",
                               "Lasso" = "#8191be")) +
  labs(title = "",
       y = "",
       x = "AUC\n(Genus level)",
       fill = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        axis.ticks.x = element_line(color = "black")
        )
p_performance_genus



# plot SH level prediction
res_pre_sh <- read.csv("table/res_prediction_SH.csv")
res_pre_sh <- res_pre_sh %>%
  filter(model %in% c("lasso", "rf")) %>%
  mutate(mechine_model  = str_replace_all(model, c( "svm" = "SVM", "rf" = "Random forest", "lasso" = "Lasso")))
res_pre_sh$cohort <- factor(res_pre_sh$cohort, levels = cohort_order)
res_pre_sh$mechine_model = factor(res_pre_sh$mechine_model, levels = rev(c("Random forest", "SVM", "Lasso")))

y_levels <- unique(res_pre_sh$cohort)
rect_data <- data.frame(y_idx = seq_along(y_levels)) %>% filter(y_idx %% 2 == 0)
p_n_feature_sh <- res_pre_sh %>% select(cohort, n_features) %>% unique() %>%
  ggplot(aes(y = cohort, x = n_features)) +
  geom_rect(data = rect_data,
            aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
  geom_col(fill = "#b2b6c1", color = "white", linewidth = 0.3) +
  labs(title = "",
       y = "",
       x = "Number of feature\n(SH level)"
       ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(color = "black", size = 8),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.ticks = element_line(color = "black")
        )
p_n_feature_sh

p_performance_sh <- res_pre_sh %>%
  ggplot(aes(y = cohort, x = auc, fill = mechine_model)) +
  geom_rect(data = rect_data,
            aes(ymin = y_idx - 0.5, ymax = y_idx + 0.5,
                xmin = -Inf, xmax = Inf),
                fill = "gray95",
                inherit.aes = FALSE) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(xmin = low, xmax = high),
                position = position_dodge(width = 0.8),
                width = 0, linewidth = 0.1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "#212121", lwd = 0.5, alpha = 0.7) +
  scale_fill_manual(values = c("SVM" = "#8cb8dc",
                               "Random forest" = "#bfdcdc",
                               "Lasso" = "#8191be")) +
  labs(title = "",
       y = "",
       x = "AUC\n(SH level)",
       fill = "Model") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        axis.ticks.x = element_line(color = "black")
        )
p_performance_sh



### part6 plot arrangement
plot1 <- p_sample_number + plot[["Observed"]] +plot[["Shannon"]] + plot[["Simpson"]] + p_adonis_robust_sh + p_adonis_robust_genus  +
    plot_layout(nrow = 1) + theme(legend.position = "right")

plot3 <- p_n_feature_sh + p_performance_sh + p_n_feature_genus + p_performance_genus +
  plot_layout(guides = "collect", nrow = 1)

final_plot <- plot1 / p_sig_feature / p_sig_feature_summary / plot3 + plot_layout(heights = c(5, 5, 1, 5))
ggsave(plot = final_plot, "plot/final_plot.pdf", width = 14, height  = 13)
