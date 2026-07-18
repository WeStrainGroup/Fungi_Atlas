# Purpose: Plot geographic distance-decay, MDS density, continent contrasts, edge sharing, heatmaps, and networks.
# Panels: Figure 1C; Figure 3F-H; Figure 4B-C; Figure S6; Figure S8A and S8C-D.

library(tidyverse)
library(phyloseq)
library(microbiome)
library(data.table)
library(vegan)
library(ggrepel)
library(ggpubr)
library(scales)
library(gridExtra)
library(paletteer)
library(igraph)
library(ggraph)
library(tidygraph)
library(ComplexUpset)
library(eulerr)
library(patchwork)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggbeeswarm)
library(broom)
library(ggridges)
library(maps)

# --- Relationship between community distance and geographical distance ---
dis <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/distance/Distance_association_rawdata_sample20.csv")
head(dis)

# compute regression metrics
model <- lm(robust_aitchison_dis ~ geographical_dis_km, data = dis)
gl <- broom::glance(model)

n_sample = gl$df.residual + 1
intercept = summary(model)$coefficients[1]
slope = summary(model)$coefficients[2]
r_squared = gl$r.squared
r_squared_adj = gl$adj.r.squared
p_lm = gl$p.value

# plot
p <- ggplot(dis, aes(x = geographical_dis_km, y = robust_aitchison_dis)) +
        geom_point(alpha = 0.1, shape = 16, size = 0.5, color = "grey70") +
        geom_smooth(method = "lm", se = F, color = "#e94f37", linewidth = 0.75) +
        # regression equation
        stat_regline_equation(
            aes(label = ..eq.label..),
            label.x = 0,
            label.y = 24.5,   # slightly below the top
            output.type = "expression") +
        # adjusted R^2 and P
        annotate("text",
                x = 0, y = 22,
                label = paste0("italic(R)[adj]^2 == ", format(r_squared_adj, digits=2),
                               " * ',' ~ italic(P) == ", format(p_lm, digits=2)),
                hjust = 0,
                parse = TRUE) +
        # Spearman rho and P
        stat_cor(method = "spearman",
              aes(label = paste0("italic(rho) == ", ..r.., " * ',' * ", " ~ italic(P) == ", ..p..)),
              label.x = 0,
              label.y = 19.5,
              parse = TRUE) +
        labs(x = "Geographical distance (km)", y = "Robust aitchison distance") +
        ylim(0,25) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", linewidth = 1))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "scatter_distance.pdf"), p, width = 3, height = 3)




# --- DB index distributions (Observed vs. Random) ---
df_DB_index <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/distance/DB_index_region.csv", row.names = 1)
head(df_DB_index)

# ensure grouping label is a factor for consistent colors
df_DB_index$Label <- factor(df_DB_index$Label)

# subset random data
random_df <- df_DB_index %>% filter(Type == "Random")
head(random_df)

random_df_stat <- random_df %>%
  filter(Type == "Random") %>%
  group_by(Label) %>%
  summarize(
    median_DB_index = median(DB_index),
    min_DB_index    = min(DB_index)
  )
random_df_stat

DB_index_median_continent <- random_df_stat$median_DB_index[random_df_stat$Label == "Continent"]
DB_index_median_country <- random_df_stat$median_DB_index[random_df_stat$Label == "Country"]
DB_index_median_region <- random_df_stat$median_DB_index[random_df_stat$Label == "Region"]

DB_index_median_continent
DB_index_median_country
DB_index_median_region

# extract observed DB_index
observed_df <- df_DB_index %>% filter(Type == "Observed")
DB_index_continent <- observed_df$DB_index[observed_df$Label == "Continent"]
DB_index_country <- observed_df$DB_index[observed_df$Label == "Country"]
DB_index_region <- observed_df$DB_index[observed_df$Label == "Region"]

DB_index_continent
DB_index_country
DB_index_region

col_continent <- "#dad7cd"
col_country <- "#a3b18a"
col_region <- "#588157"

col <- c("Continent" = col_continent,
         "Country" = col_country,
         "Region" = col_region)

# plot
his <- ggplot(random_df, aes(x = DB_index, y = Label, fill = Label)) +
  geom_density_ridges(
    alpha = 0.5,
    scale = 1.5,        # ridge height
    color = "black",
    rel_min_height = 0
  ) +
  scale_fill_manual(values = col) +
  geom_segment(aes(x = DB_index_continent, xend = DB_index_continent, y = 1, yend = 2),
             color = col_continent, linewidth = 0.5) +
  geom_segment(aes(x = DB_index_country, xend = DB_index_country, y = 2, yend = 3),
             color = col_country, linewidth = 0.5) +
  geom_segment(aes(x = DB_index_region, xend = DB_index_region, y = 3, yend = 4.45),
             color = col_region, linewidth = 0.5) +
  annotate("text", x = 160, y = 4.1,
           label = paste0("level: region",
                          "\ntrue value: ", round(DB_index_region, 1)
                           #"\nrandom median: ", round(DB_index_mean_region, 1)
                           ),
           size = 3, hjust = 0) +
  annotate("text", x = 160, y = 2.65,
           label = paste0("level: country",
                          "\ntrue value: ", round(DB_index_country, 1)
                           #"\nrandom median: ", round(DB_index_mean_country, 1)
                           ),
           size = 3, hjust = 0) +
  annotate("text", x = 160, y = 1.65,
           label = paste0("level: continent",
                          "\ntrue value: ", round(DB_index_continent, 1)
                           #"\nrandom median: ", round(DB_index_mean_continent, 1)
                           ),
           size = 3, hjust = 0) +
  labs(
    x = "Davies-Bouldin's index",
    y = "Iterations",
    fill = NULL,
    color = NULL
  ) +
  xlim(0, 250) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.key.size = unit(0.15, "in"),
            legend.position = "none")


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "his_DB_index.pdf"), his, width = 3, height = 3)


# --- PCoA visualization ---
meta <- read.csv("ID_meta_final.csv", header = T, row.names = 1)
res <- read.csv("mds_reduction/res_MDS.csv", row.names = 1)
res <- res %>% left_join(meta, by = "ID")
explained <- read.csv("mds_reduction/res_MDS_explained.csv", row.names = 1)

# plot by continent
Continent_list <- c( "Asia", "North America","Europe",  "Africa")
colors <- c("white", "grey", "#5A4D9F", "#5C6CA0", "#5F8BA2",
            "#62AAA3", "#82C383", "#B1DA51", "#DFF020", "#EFE802",
            "#C0A408", "#92610F", "#641E16", "#632929")
p <- list()
for(i in 1:4){
  p[[i]] <- res %>%
    filter(Continent == Continent_list[i]) %>%
    ggplot(aes(x = V1, y =  V2)) +
    geom_point(size = 0.01, color = "#606060")+
    geom_density_2d_filled(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # x = 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # y = 0
    scale_fill_manual(values = colors) +
    xlim(-2000, 2000)+
    ylim(-400, 400)+
    coord_fixed(ratio = 40/8) +
    labs(x = "", y = "")+
    annotate("text", x = -2000, y = 400, label = Continent_list[i], hjust = 0, vjust = 1) +
    theme_minimal()+
    theme(
      legend.position = "non",
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
    )
}
plot_continent <- p[[1]] + p[[2]] + p[[3]] + p[[4]] + plot_layout(nrow = 1)
plot_continent


# plot by country
Country_summary <- table(res$Country) %>% data.frame() %>% arrange(-Freq)
Country_list <- Country_summary$Var1
length(Country_list)
p <- list()
for(i in 1:24){
  p[[i]] <- res %>%
    filter(Country == Country_list[i]) %>%
    ggplot(aes(x = V1, y =  V2)) +
    geom_point(size = 0.01, color = "#606060")+
    geom_density_2d_filled(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal()+
    xlim(-2000, 2000)+
    ylim(-400, 400)+
    coord_fixed(ratio = 40/8) +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "")+
    annotate("text", x = -2000, y = 400, label = Country_list[i], hjust = 0, vjust = 1) +
    theme(
      legend.position = "non",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
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
    )
}
plot_country <- p[[1]] + p[[2]] + p[[3]] + p[[4]] + p[[5]] + p[[6]] + p[[7]] +p[[8]] + p[[9]] + p[[10]] +
  p[[11]] + p[[12]] + p[[13]] + p[[14]] + p[[15]] + p[[16]] + p[[17]] + p[[18]] + p[[19]] + p[[20]] +
  p[[21]] + p[[22]] + p[[23]] + p[[24]] + plot_layout(ncol = 4)
plot_country




# --- Geographic distribution of sample origins ---
setwd("data_for_analysis/")

sample_distri <- read.csv("sample_info/plotdata_sample_distribution.csv")
country_order <- sample_distri %>% select(Country, Continent) %>% unique() %>% arrange(Continent, Country)
sample_distri$Country <- factor(sample_distri$Country, levels = country_order$Country)
colors <- c(
  "Burkina Faso"             = "#6c838b", #Africa
  "Central African Republic" = "#9a6225",
  "Madagascar"               = "#e49edd",
  "China"                    = "#cd5c5c", #Asia
  "India"                    = "#fbbdd4",
  "Israel"                   = "#29b95c",
  "Japan"                    = "#4463d8",
  "Malaysia"                 = "#00c5b6",
  "Republic of Korea"        = "#f6c6ac",
  "Thailand"                 = "#975ccb",
  "United Arab Emirates"     = "#34495e",
  "Austria"                  = "#901eb4", #Europe
  "Belgium"                  = "#d9f2d0",
  "Czechia"                  = "#f1a984",
  "Denmark"                  = "#aaffc3",
  "France"                   = "#fbe5d5",
  "Germany"                  = "#ffe019",
  "Ireland"                  = "#ff69b4",
  "Italy"                    = "#688620",
  "Lithuania"                = "#1399b2",
  "Netherlands"              = "#7f6000",
  "Norway"                   = "#60cbf3",
  "Spain"                    = "#844b44",
  "Canada"                   = "#bfee46", #North America
  "United States of America" = "#ff8c69",
  "Mexico"                   = "#017a73"
)

# create world background
world_background <- ggplot() +
  geom_polygon(data = map_data("world"),aes(x=long, y=lat, group=group),fill="darkgrey")+
  theme_bw()+
  scale_x_continuous(expand = expansion(add=c(0,0)))

# add sample to the world background
sample_map <- world_background +
  geom_point(data=sample_distri %>% filter(!is.na(Country)),
             aes(x=Longitude,y=Latitude, fill=Country),
             shape = 21, size=2, color= "black", alpha = 0.85) +
  #geom_text(data=sample_distri, aes(x=Longitude,y=Latitude, label = Country), vjust = 0, color = "black", size = 2)+
  scale_fill_manual(values = colors)+
  ylim(-66, 90)+
  labs(x = "Longitude", y = "Latitude", fill = "") +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.position = "left",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.ticks.length = unit(1, "mm"),
    legend.text = element_text(lineheight = 0.4, size = 6),
    legend.key.height = unit(0.3, "cm")
  )
sample_map
ggsave(plot = sample_map, "plot/Global_sample_map.pdf", width = 7, height = 3)




# --- Upset plot for unique/shared genera across continents ---
gen_overlap <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/compare/overlap_gen_0.01_filterdp.csv",
                        row.names = 1)
head(gen_overlap)

# remove rows containing Other|Unknown
gen_overlap <- gen_overlap[!grepl("Other|Unknown", rownames(gen_overlap)), ]

# build input dataframe
upset_df <- gen_overlap[, c("Africa", "Asia", "Europe", "NorthA")]

# base upset plot
upset <- upset(upset_df,
               intersect = c("Africa", "Asia", "Europe", "NorthA"),
               name = "Genus Overlap Across Continents (Prevalence > 1%)",

               base_annotations = list('Num of genus' = intersection_size(text = list(size = 3))),

              # add set size labels
              set_sizes = (upset_set_size() +
                           geom_bar(width = 0.8) +   # adjust bar width
                           geom_text(aes(label = ..count..),
                                     stat = "count",
                                     hjust = -0.1,
                                     size = 3,
                                     color = "white") +
                                     ylab("Num of genus")),

              width_ratio = 0.25, # control ratio of panels
              sort_intersections_by = "cardinality",

)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "upset_gen.pdf"), upset, width = 6, height = 3)



# --- Visualizing continent-level differential analysis results ---
setwd("data_for_analysis/")
continent_list <- c("Asia", "North America", "Europe", "Africa")

mean_abun <- read.csv("abundance_prevalence/genus_p0.01_log10relative_abundance_mean_countinent.csv", row.names = 1) %>% rename(`North America` = North.America)
metadata <- read.csv("ID_meta_final.csv", sep = ",", header = T, row.names = 1) #%>% rename(Sample_ID = ID)
prevalence <- read.csv("abundance_prevalence/genus_p0.01_prevelance.csv", row.names = 1)
res <- read.csv("differential_analysis_res/plotdata_continent_diff_genus.csv", row.names = 1)

res_pre <- res %>%
  left_join(prevalence %>% rename_with(~c("Taxa", "Variable_pre", "Variable")), by = c("Taxa", "Variable")) %>%
  left_join(prevalence %>% rename_with(~c("Taxa", "Reference_pre", "Reference_variable")), by = c("Taxa", "Reference_variable")) %>%
  filter(Significance1 == "Sig") %>%
  select(Variable, Reference_variable, Taxa, Estimate, adjust.p.value, Variable_pre, Reference_pre, Significance2) %>%
  data.frame()  %>%
  filter(Significance2 == "Higher") # select higher in variable

# file
res2 <- res %>% filter(Significance2 == "Higher") %>%
  mutate(Label = paste0(Variable, " vs. ", Reference_variable)) %>%
  mutate(taxa = str_replace_all(Taxa, "g__", "")) %>%
  select(Label, taxa, p.value, adjust.p.value) %>%
  mutate(significance = case_when(is.na(adjust.p.value) ~ 0, T ~ 1))
res3 <- res2 %>% select(Label, taxa, significance) %>%
  pivot_wider(names_from = "taxa", values_from = "significance") %>%
  pivot_longer(cols = unique(res2$taxa), names_to = "taxa", values_to = "significance") %>%
  mutate(significance = case_when(is.na(significance) | significance == 0 ~ "non-Sig", T ~ "Sig")) %>%
  mutate(temp = Label) %>% separate(col = temp, into = c("Continent", "Reference"), sep = " vs. ") %>%
  mutate(label2 = case_when(significance == "Sig" ~ Continent, T ~ "non-Sig"))

mean_abun_long <- mean_abun %>% rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = all_of(continent_list), names_to = "Continent", values_to = "Relative_abundance")
abun1 <- mean_abun_long %>% mutate(taxa = str_replace_all(Taxa, "g__", "")) %>% filter(taxa %in% res2$taxa)
abun2 <- abun1 %>% mutate(abundance = Relative_abundance) %>%
  select(taxa, abundance) %>%
  group_by(taxa) %>%
  summarise(abundance_mean = mean(abundance, na.rm = T), .groups = "drop") %>%
  ungroup %>%
  arrange(abundance_mean)
taxa_list <- abun2$taxa

pre1 <- prevalence %>% mutate(taxa = str_replace_all(Genus, "g__", "")) %>% filter(taxa %in% res2$taxa)
continent_list <- c("Asia",  "North America", "Europe", "Africa")
color_continent <- c("Asia" =          "#f09ba0",
                     "North America" = "#1772f6",
                     "Europe" =        "#4d0086",
                     "Africa" =        "#edb11f")


# data for left figure
compare_order <- read.table("differential_analysis_res/continent_compare_order.txt", sep = "\t", header = T, check.names = F)
compare_order_long <- compare_order %>%
  pivot_longer(cols = c("Asia", "North America", "Europe", "Africa"), names_to = "continent", values_to = "label") %>%
  filter(!is.na(label)) %>%
  mutate(label = as.factor(label), order = as.factor(order)) %>%
  arrange(order)
compare_order_long$order <- factor(compare_order_long$order, levels = rev(unique(compare_order_long$order)))
compare_order_long$continent <- factor(compare_order_long$continent, levels = c("Asia",  "North America", "Europe", "Africa"))
compare_order_long <- compare_order_long %>% mutate(label2 = case_when(label == "1" ~ continent, T ~ "Reference continent"))

# data for right figure
compare_summary <- table(res3$Label[res3$significance == "Sig"]) %>%
  data.frame() %>% rename(Label = Var1)
compare_summary$Continent <- c("Africa", "Africa","Africa",
                               "Asia", "Asia","Asia",
                               "Europe","Europe","Europe",
                               "North America","North America","North America")
compare_list <- c(
  "Asia vs. North America", "Asia vs. Europe", "Asia vs. Africa",
  "North America vs. Asia", "North America vs. Europe", "North America vs. Africa",
  "Europe vs. Asia", "Europe vs. North America", "Europe vs. Africa",
  "Africa vs. Asia",  "Africa vs. North America", "Africa vs. Europe"
  )
compare_summary$Label <- factor(compare_summary$Label, levels = rev(compare_list))

# left figure
p_compare_order <-
  ggplot(data = compare_order_long, aes(x = continent, y = order)) +
  geom_point(aes(fill = label2), size = 4, shape = 21, alpha = 1)+
  scale_fill_manual(values = c(color_continent, "Reference continent" = "white"))+
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

# mid figure
res3$Label <- factor(res3$Label, levels = rev(compare_list))
#res3$taxa <- factor(res3$taxa, levels = taxa_list)
p_compare <-
  ggplot(res3, aes(x = taxa, y = Label, fill = label2, alpha = 0.85)) +
  geom_tile(color = "white", alpha = 1, linewidth = 0.5) +
  scale_fill_manual(values = c(color_continent, "non-Sig" = "#ebebeb")) +
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
    #axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks.length.y = unit(0, "mm"),
  )
p_compare
# right figure
p_compare_summary <- compare_summary %>%
  ggplot(aes(x = Freq, y = Label, fill = Continent)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(x = Freq, label = Freq),
            hjust = 0.5,
            size = 3,
            nudge_x = 1,
            color = "black") +
  xlim(0, 25)+
  labs(y = "", x = "Number of specific taxa") +
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

# figure: left + mid + right
p_total <- p_compare_order + p_compare + p_compare_summary + plot_layout(nrow = 1, widths = c(1,  10, 2))


### supplementary figure
# figure: relative abundance
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

# figure: prevalence
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
p_prev_abun <- p_abun + p_pre[[1]] + p_pre[[2]] + p_pre[[3]] + p_pre[[4]] + plot_layout(nrow = 1, widths = c(10, 2, 2, 2, 2))






# --- Unique and shared significant edges across continents ---
edge_overlap <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/compare/overlap_edge_0.01_filterdp.csv",
                        row.names = 1)
head(edge_overlap)
nrow(edge_overlap)

# prepare upset data
upset_df <- edge_overlap[, c("Africa_presence", "Asia_presence", "Europe_presence", "NorthA_presence", "direction")]

# remove edges present in only one region (direction-specific uniques)
upset_df <- upset_df %>%
  filter(Africa_presence + Asia_presence + Europe_presence + NorthA_presence > 1)

# rename columns
upset_df <- upset_df %>%
  rename(
    Africa = Africa_presence,
    Asia = Asia_presence,
    Europe = Europe_presence,
    NorthA = NorthA_presence
  )
head(upset_df)

# define colors
col <- c(
  positive = scales::alpha("#f87060", 0.9),
  negative = scales::alpha("#102542", 0.9),
  mixed = scales::alpha("#cdd7d6", 0.9)
)


upset <- upset(data = upset_df,
               intersect = c("Africa", "Asia", "Europe", "NorthA"),
               name = "Shared edges (p-value ≤ 0.001) across continents",
               width_ratio = 0.25,
               sort_intersections_by = "cardinality",

               base_annotations = list('Num of edge' = (intersection_size(mapping = aes(fill = direction),
                                                                text = list(size = 3)) +
                                              scale_fill_manual(values = col, name = "Direction") +
                                              theme(legend.position = c(0.75, 1),
                                                    legend.justification = c(0, 1)))), # align top-left


              # add set size labels
              set_sizes = (upset_set_size() +
                           geom_bar(width = 0.8) +   # adjust bar width
                           geom_text(aes(label = ..count..),
                                     stat = "count",
                                     hjust = -0.1,
                                     size = 3,
                                     color = "white") +

                  ylab("Num of edge"))
)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "upset_edge.pdf"), upset, width = 6, height = 3)




# --- Pie chart of edge directions ---
# count by direction
direction_count <- upset_df %>%
  count(direction) %>%
  mutate(percent = n / sum(n),
         label = paste0(round(percent * 100), "%"))

# draw pie chart
pie1 <- ggplot(direction_count, aes(x = "", y = n, fill = direction)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 8, color = "black") +
  scale_fill_manual(values = col,
                    name = NULL,
                    guide = "none") +
  theme_void() +
  theme(legend.position = "none")

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "pie_edge_direction.pdf"), pie1, width = 3, height = 3, bg = "transparent")





# --- Bar plot of counts: significant edges and shared edges by region ---
edge_overlap <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/compare/overlap_edge_0.01_filterdp.csv",
                        row.names = 1)


# total significant edges per region
region_edge_counts <- edge_overlap %>%
  summarise(
    Africa = sum(Africa_presence, na.rm = TRUE),
    Asia = sum(Asia_presence, na.rm = TRUE),
    Europe = sum(Europe_presence, na.rm = TRUE),
    NorthA = sum(NorthA_presence, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "Region", values_to = "EdgeCount")

region_edge_counts


# edges appearing in at least two regions
shared_edges <- edge_overlap[edge_overlap$all_count >= 2, ]
nrow(edge_overlap) # 17474
nrow(shared_edges) # 1922

region_shared_edge_counts <- shared_edges %>%
  summarise(
    Africa = sum(Africa_presence, na.rm = TRUE),
    Asia = sum(Asia_presence, na.rm = TRUE),
    Europe = sum(Europe_presence, na.rm = TRUE),
    NorthA = sum(NorthA_presence, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "Region", values_to = "SharedEdgeCount")

region_shared_edge_counts



# merge total and shared counts
region_edge_merged <- merge(region_edge_counts, region_shared_edge_counts, by = "Region")
region_edge_merged
region_edge_merged$NotSharedEdgeCount <- region_edge_merged$EdgeCount - region_edge_merged$SharedEdgeCount
region_edge_merged

# reshape to long format and compute proportions of unique vs. shared per region
df_long <- region_edge_merged %>%
  pivot_longer(
    cols = c(SharedEdgeCount, NotSharedEdgeCount),
    names_to = "EdgeType",
    values_to = "Count") %>%
  group_by(Region) %>%
  mutate(Proportion = Count / sum(Count))

df_long

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/fastspar"
write.csv(df_long, file.path(path, "edge_uniquness_plot.csv"))

# ggplot
p <- ggplot(df_long, aes(x = Region, y = Proportion, fill = EdgeType)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
  # count labels (raw numbers)
  geom_text(
    aes(label = Count),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3.5,
    angle = 90  # rotate labels 90 degrees
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "SharedEdgeCount" = "#264653",
      "NotSharedEdgeCount" = "#2a9d8f"
    ),
    labels = c(
      "SharedEdgeCount" = "Shared",
      "NotSharedEdgeCount" = "Unique"
    ),
    name = NULL
  ) +
  labs(y = "Proportion of edges with p ≤ 0.001") +
  theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.15, "in"),
    axis.title.x = element_blank()
  )



path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_shared_edge.pdf"), p, width = 2.25, height = 4, bg = "transparent")



# --- Heatmap of edges significant in all four regions ---
# split edge into Node1 and Node2
head(edge_overlap)
edge_overlap2 <- edge_overlap %>%
  filter(all_count == 4) %>%
  select(edge, Africa_cor, Asia_cor, Europe_cor, NorthA_cor, pos_count) %>%
  rename_with(~sub("_cor$", "", .x), ends_with("_cor")) %>%
  column_to_rownames("edge")

table(edge_overlap2$pos_count)

mat <- as.matrix(edge_overlap2[, c("Africa", "Asia", "Europe", "NorthA")])
mat[1:3, 1:3]

# cluster edges within pos_count groups and order
# initialize container for final order
final_order <- c()

# get group names
groups <- sort(unique(edge_overlap2$pos_count))

for (grp in groups) {
  idx <- which(edge_overlap2$pos_count == grp)
  sub_mat <- mat[idx, , drop = FALSE]

  if (nrow(sub_mat) >= 2) {
    # compute Euclidean distance
    dist_sub <- dist(sub_mat)
    # hierarchical clustering
    hc_sub <- hclust(dist_sub)
    # extract order
    ordered_edges <- rownames(sub_mat)[hc_sub$order]
  } else {
    # only one row: keep as is
    ordered_edges <- rownames(sub_mat)
  }

  final_order <- c(final_order, ordered_edges)
}

ordered_edges <- final_order

# cluster regions (columns)
d_col <- dist(t(mat))
clust_col <- hclust(d_col)
ordered_regions <- clust_col$labels[clust_col$order]


# convert to long format
edge_long <- edge_overlap2 %>%
  rownames_to_column("edge") %>%
  separate(edge, into = c("Node1", "Node2"), sep = "-", remove = FALSE) %>%
  pivot_longer(
    cols = -c(edge, Node1, Node2, pos_count),
    names_to = "region",
    values_to = "correlation") %>%
  mutate(
    edge = factor(edge, levels = ordered_edges),  # reorder by clustering
    region = factor(region, levels = ordered_regions)
  )



# heatmap
heat <- ggplot(edge_long, aes(y = region, x = edge, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0,
    name = "correlation"
  ) +
  labs(
    y = "",
    x = "Edges that consistently significant correlated across all continents"
  ) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black", linewidth = 0.5)+
  geom_vline(xintercept = 19.5, linetype = "dashed", color = "black", linewidth = 0.5)+
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #axis.text.x = element_text(face = "italic"),
    legend.key.size = unit(0.15, "in"),
  )



path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "heat_shared_edge.pdf"), heat, width = 8, height = 3, bg = "transparent")




# --- Genus-level taxa co-occurrence network ---
location <- "Europe"
location <- "Africa"
location <- "NorthA"
location <- "Asia"
location <- "All"

file <- file.path("/data/wangxinyu/ITS_Public/Project/Analysis/fastspar",
 paste0(location, "_gen_0.01"),
 paste0(location, "_gen_0.01_edge.csv"))
edge <- read.csv(file, row.names = 1)
head(edge)

# filter edges
top <- 50
threshold_p <- 0.001

edge_filter <- edge %>%
  # exclude "Other" and "Unknown" and filter by p-value
  filter(
    !grepl("Other|Unknown", node1) &
    !grepl("Other|Unknown", node2),
    pval <= threshold_p
  ) %>%
  # take top by |cor|
  arrange(desc(abs(cor))) %>%
  slice_head(n = top)

edge_filter


# build network object
net <- graph_from_data_frame(edge_filter, directed = FALSE) %>%
  as_tbl_graph() %>%
  activate(edges) %>%
  mutate(
    strength = abs(cor),
    direction_flag = cor > 0  # TRUE = positive, FALSE = negative
  ) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree())

# check maxima
max_degree <- net %>%
  activate(nodes) %>%
  as_tibble() %>%
  pull(degree) %>%
  max(na.rm = TRUE)
max_degree #Europe 10, Africa 12, NorthA 11, Asia 10, All 15


max_strength <- net %>%
  activate(edges) %>%
  as_tibble() %>%
  pull(strength) %>%
  max(na.rm = TRUE)
max_strength #Europe 0.37, Africa, 0.64, NorthA, 0.48, Asia 0.48, All 0.5279

min_strength <- net %>%
  activate(edges) %>%
  as_tibble() %>%
  pull(strength) %>%
  min(na.rm = TRUE)
min_strength #Europe 0.2096, Africa, 0.2924, NorthA, 0.1604, Asia 0.1181, All 0.1276


# add phylum metadata
meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/tax_gen_filterdp.csv", row.names = 1)
meta <- meta %>% select(Phylum)

# remove prefixes
rownames(meta) <- sub("^[a-z]{1,2}__", "", rownames(meta))
meta$Phylum <- sub("^[a-z]{1,2}__", "",meta$Phylum )
head(meta)

net <- net %>%
  activate(nodes) %>%
  mutate(taxon = name) %>%  # rename 'name' to 'taxon' for join
  left_join(
    meta %>%
      rownames_to_column("taxon"),   # move rownames to column for matching
    by = "taxon"
  )


# custom edge color mapping: TRUE=red, FALSE=blue
edge_colors <- c("TRUE" = "#DE6449", "FALSE" = "#4F86C6")
phylum_colors <- c(Ascomycota = "#333030",
                   Basidiomycota = "#9055A2",
                   Mucoromycota = "#4f953b",
                   Chytridiomycota = "#c1121f",
                   Mortierellomycota = "#457b9d")


# plot without legend
p_net <- ggraph(net, layout = "circle") +
  geom_edge_link(aes(width = strength, color = as.character(direction_flag)), alpha = 0.5) +
  geom_node_point(aes(size = degree, color = Phylum), shape = 16, alpha = 0.95) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +

  # manual node color mapping
  scale_color_manual(
    values = phylum_colors,
    limits = names(phylum_colors),  # ensure all phyla appear in legend
    name = NULL
  ) +
  # manual edge color mapping
  scale_edge_color_manual(
    values = edge_colors,
    breaks = c("TRUE", "FALSE"),
    labels = c("positive", "negative"),
    name = NULL
  ) +
  # manual edge width mapping
  scale_edge_width(name = "correlation", limits = c(0.1, 0.65)) +
  # manual node size mapping
  scale_size_continuous(name = "degree", limits = c(1, 12)) +
  # legend node size override
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "none")


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, paste0("net_", location, ".pdf")), p_net, width = 4, height = 4)


# export legend only
# Europe top1000
legend_only <- ggraph(net, layout = "circle") +
  geom_edge_link(aes(width = strength, color = as.character(direction_flag)), alpha = 0.5) +
  geom_node_point(aes(size = degree, color = Phylum), shape = 16, alpha = 0.95) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +

  # legend scales
  scale_color_manual(
    values = phylum_colors,
    limits = names(phylum_colors),
    name = "Phylum",
    drop = FALSE # show all colors
  ) +
  scale_edge_color_manual(
    values = edge_colors,
    breaks = c("TRUE", "FALSE"),
    labels = c("positive", "negative"),
    name = "Correlation direction"
  ) +
  scale_edge_width(name = "Correlation strength", limits = c(0.1, 0.65)) +
  scale_size_continuous(name = "Degree", limits = c(1, 12)) +

  # force each legend to a single row; order controls sequence
  guides(
    color = guide_legend(nrow = 1, order = 1, override.aes = list(size = 4)),
    edge_color = guide_legend(nrow = 1, order = 2),
    edge_width = guide_legend(nrow = 1, order = 3),
    size = guide_legend(nrow = 1, order = 4)
  ) +

  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",  # stack legend blocks vertically, each on one row
    legend.box.just = "center",
    legend.title = element_text(size = 12),     # legend title size
    legend.text = element_text(size = 12),      # legend text size
  plot.margin = margin(t = 5, r = 5, b = 25, l = 5)  # increase bottom margin
  )

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "net_legend.pdf"), legend_only, width = 16, height = 4)


# standalone plot for 'All' with legend
p_net2 <- ggraph(net, layout = "circle") +
  geom_edge_link(aes(width = strength, color = as.character(direction_flag)), alpha = 0.5) +
  geom_node_point(aes(size = degree, color = Phylum), shape = 16, alpha = 0.95) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = phylum_colors, name = NULL) +
  scale_edge_color_manual(values = edge_colors, breaks = c("TRUE", "FALSE"),
                          labels = c("positive", "negative"),
                          name = NULL,
                          guide = guide_legend(override.aes = list(linewidth = 2))) + # name=NULL removes legend title
  scale_edge_width(name = "correlation") +
  scale_size_continuous(name = "degree") +
  guides(color = guide_legend(override.aes = list(size = 4))) + # thicker line in legend
  theme_void()

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "net_All.pdf"), p_net2, width = 5, height = 4)
