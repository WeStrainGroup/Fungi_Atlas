# Purpose: Analyze distance-decay, MDS, cluster strength, continent-enriched taxa, and shared correlation edges.
# Panels: Figure 3F-H; Figure 4A-C; Figure S6; Figure S8A-D.

library(phyloseq)
library(microbiome)
library(data.table)
library(Biostrings)
library(tidyverse)
library(vegan)
library(effectsize)
library(geosphere)
library(bigmds)
library(broom)
library(future)
library(future.apply)
library(clusterSim)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(iNEXT)
library(parallel)
library(rstatix)

# Read a large OTU table and restore its first column as row names.
read_otu <- function(filepath) {
  otu <- fread(filepath, header = TRUE)
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[1]]
  otu <- otu[, -1]
  return(otu)
}

# --- Association between fungal community distance and geographic distance ---
etwd("data_for_analysis/")
otu <-  readRDS("otu_gen_0.001_robait_filterdp_dist.rds")
otu_dis <- as.matrix(otu)
metadata <- read.csv("ID_meta_final.csv", row.names = 1)
geo_dis <- read.table("distance_association/geographic_distance.tsv", sep = "\t", header = T)
pro_list <- metadata$Region %>% unique()
pro_list1 <- pro_list[!str_detect(pro_list, "Unknown")] %>% setdiff("missing")
id_list1 <- metadata$Sample_ID[which(metadata$Region %in% pro_list1)] %>%
  intersect(row.names(otu_dis))

# Distance calculation
meta <- read.csv("ID_meta_final.csv")
df1 <- meta %>% select(Region, Longitude, Latitude) %>%
  filter(!Region == "NA_NA") %>%
  group_by(Region) %>%
  slice(1) %>%
  ungroup() %>%
  column_to_rownames(var = "Region") %>%
  as.matrix()

distance_matrix <- distm(df1, fun=distGeo)
distance_matrix_km <- round(distance_matrix / 1000, 2) # m to km
rownames(distance_matrix_km) <- rownames(df1)
colnames(distance_matrix_km) <- rownames(df1)

distance <- distance_matrix_km %>%
  as.data.frame() %>%
  rownames_to_column(var = "Country1") %>%
  pivot_longer(cols = -1, names_to = "Country2", values_to = "Distance_km") %>%
  mutate(Distance_km = round(Distance_km, 2)) %>%
  data.frame()

write.table(distance, "distance_association/geographic_distance.tsv", sep = "\t", quote = F, row.names = F)


# Associate distances
a <- metadata %>% filter(!str_detect(Region, "Unknown") & !Country == "NA")
b <- table(a$Country)%>% data.frame() %>% arrange(-Freq)
cut_off=20 # Randomly sample 'cut_off' IDs per country
metadata_select_1 <- a %>% filter(Country %in% b$Var1[b$Freq >= cut_off]) %>%
  group_by(Country) %>%
  slice_sample(n = cut_off) %>%
  ungroup()
metadata_select_2 <- a %>% filter(Country %in% b$Var1[b$Freq < cut_off])
metadata_select <- rbind(metadata_select_1, metadata_select_2)

id_list_select <- metadata_select$ID[which(!str_detect(metadata_select$Region, "unknown"))] %>%
  intersect(row.names(otu_dis))
dis_data_select <- otu_dis[id_list_select, id_list_select] %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(cols = all_of(id_list_select), names_to = "Sample2", values_to = "Archison_distance") %>%
  mutate(Country1 = metadata$Region[match(Sample1, metadata$ID)],
         Country2 = metadata$Region[match(Sample2, metadata$ID)]) %>%
  left_join(geo_dis, by = c("Country1", "Country2"))

dis_data0 <- dis_data_select %>% data.frame() %>%
  rename(location1 = Country1, location2 = Country2)
spearman_cor <- cor.test(dis_data0$Archison_distance, dis_data0$Distance_km, method = "spearman")

write.csv(dis_data0, "distance_association/Distance_association_rawdata_sample20.csv", row.names = F)
write_rds(spearman_cor, "distance_association/spearman_correlation.rds")





# --- Geography-related PCoA ---
setwd("data_for_analysis/")
otu_gen_dist <- readRDS("otu_gen_0.001_robait_filterdp_dist.rds")
meta <- read.csv("ID_meta_final.csv", header = T, row.names = 1)

arc_dis <- as.matrix(otu_gen_dist)
set.seed(123)
mds <- divide_conquer_mds(x = arc_dis, l = 300, c_points = 5 * 10, r = 10, n_cores = 16)
explained <- data.frame(Variable = paste0("V", 1:10), Ratio = mds$eigen / sum(mds$eigen))

res <- data.frame(mds$points)
rownames(res) <- rownames(arc_dis)
res <- res %>% rename_with(~paste0("V", 1:10)) %>% rownames_to_column(var = "ID")
explained

write.csv(res, "mds_reduction/res_MDS.csv")
write.csv(explained, "mds_reduction/res_MDS_explained.csv")




# --- Evaluate clustering tightness across regions ---
# Prepare dataset
MDS <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/distance/res_MDS.csv", row.names = 1)
head(MDS)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
meta_filter <- meta %>% select("ID", "Continent", "Country", "Region")
head(meta_filter)
unique(meta_filter$Continent)
unique(meta_filter$Country)
unique(meta_filter$Region)


MDS_meta_filter <- merge(MDS, meta_filter, by = "ID", all.x = T)
head(MDS_meta)
table(MDS_meta$Region, useNA = "ifany")


# Start calculation
set.seed(6311)
n_iter <- 100000
label_cols <- c("Continent", "Country", "Region")

all_db <- data.frame(Label=character(), Iteration=integer(), DB_index=numeric(), Type=character(), stringsAsFactors=FALSE)

for (label in label_cols) {
  cat("Processing", label, "\n")

  # Filter NA and 'NA_NA'
  df_tmp <- MDS_meta_filter %>%
  filter(
    !is.na(.data[[label]]),           # Exclude NA
    .data[[label]] != "NA_NA"         # Exclude "NA_NA"
  )

  # Feature matrix
  X_tmp <- df_tmp %>% select(-ID, -Continent, -Country, -Region) %>% as.matrix()

  # Encode labels as integers
  cl_tmp <- as.numeric(factor(df_tmp[[label]]))

  # Observed DB index
  db_true <- index.DB(X_tmp, cl_tmp, centrotypes="centroids")$DB

  # Randomized DB indices
  db_random <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    cl_perm <- sample(cl_tmp)
    db_random[i] <- index.DB(X_tmp, cl_perm, centrotypes="centroids")$DB
  }

  # Combine results
  all_db <- rbind(
    all_db,
    data.frame(Label=label, Iteration=0, DB_index=db_true, Type="Observed"),
    data.frame(Label=label, Iteration=1:n_iter, DB_index=db_random, Type="Random")
  )
}


head(all_db)
nrow(all_db)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/distance"
write.csv(all_db, file.path(path, "DB_index_region.csv"))









# --- Split data by four continents ---
otu <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/otu_ASV_filterdp.csv")
otu[1:3, 1:3]

taxa <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/tax_ASV_filterdp.csv", row.names = 1)
head(taxa)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
unique(meta$Continent)
table(meta$Continent)

rownames(meta) <- meta$ID
head(meta)

ps <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = FALSE),
                     tax_table(as.matrix(taxa)),
                     sample_data(meta))
ps

ps_Asia <- subset_samples(ps, Continent == "Asia")
ps_Africa <- subset_samples(ps, Continent == "Africa")
ps_Europe <- subset_samples(ps, Continent == "Europe")
ps_NorthA <- subset_samples(ps, Continent == "North America")

# Remove taxa absent in each subset
ps_Asia <- prune_taxa(taxa_sums(ps_Asia) > 0, ps_Asia)
ps_Africa <- prune_taxa(taxa_sums(ps_Africa) > 0, ps_Africa)
ps_Europe <- prune_taxa(taxa_sums(ps_Europe) > 0, ps_Europe)
ps_NorthA <- prune_taxa(taxa_sums(ps_NorthA) > 0, ps_NorthA)

ps_Asia
ps_Africa
ps_Europe
ps_NorthA

ps_Asia_gen <- aggregate_taxa(ps_Asia, level = "Genus")
ps_Africa_gen <- aggregate_taxa(ps_Africa, level = "Genus")
ps_Europe_gen <- aggregate_taxa(ps_Europe, level = "Genus")
ps_NorthA_gen <- aggregate_taxa(ps_NorthA, level = "Genus")

ps_Asia_gen
ps_Africa_gen
ps_Europe_gen
ps_NorthA_gen

ps_Asia_gen_0.01 <- aggregate_rare(ps_Asia, level = "Genus", detection = 0, prevalence = 0.01)
ps_Africa_gen_0.01 <- aggregate_rare(ps_Africa, level = "Genus", detection = 0, prevalence = 0.01)
ps_Europe_gen_0.01 <- aggregate_rare(ps_Europe, level = "Genus", detection = 0, prevalence = 0.01)
ps_NorthA_gen_0.01 <- aggregate_rare(ps_NorthA, level = "Genus", detection = 0, prevalence = 0.01)

ps_Asia_gen_0.01
ps_Africa_gen_0.01
ps_Europe_gen_0.01
ps_NorthA_gen_0.01

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific"
write.csv(otu_table(ps_Asia), file = file.path(path, "Asia_otu_ASV_filterdp.csv"))
write.csv(tax_table(ps_Asia), file = file.path(path, "Asia_tax_ASV_filterdp.csv"))
write.csv(otu_table(ps_Africa), file = file.path(path, "Africa_otu_ASV_filterdp.csv"))
write.csv(tax_table(ps_Africa), file = file.path(path, "Africa_tax_ASV_filterdp.csv"))
write.csv(otu_table(ps_Europe), file = file.path(path, "Europe_otu_ASV_filterdp.csv"))
write.csv(tax_table(ps_Europe), file = file.path(path, "Europe_tax_ASV_filterdp.csv"))
write.csv(otu_table(ps_NorthA), file = file.path(path, "NorthA_otu_ASV_filterdp.csv"))
write.csv(tax_table(ps_NorthA), file = file.path(path, "NorthA_tax_ASV_filterdp.csv"))

write.csv(otu_table(ps_Asia_gen), file = file.path(path, "Asia_otu_gen_p0_filterdp.csv"))
write.csv(tax_table(ps_Asia_gen), file = file.path(path, "Asia_tax_gen_p0_filterdp.csv"))
write.csv(otu_table(ps_Africa_gen), file = file.path(path, "Africa_otu_gen_p0_filterdp.csv"))
write.csv(tax_table(ps_Africa_gen), file = file.path(path, "Africa_tax_gen_p0_filterdp.csv"))
write.csv(otu_table(ps_Europe_gen), file = file.path(path, "Europe_otu_gen_p0_filterdp.csv"))
write.csv(tax_table(ps_Europe_gen), file = file.path(path, "Europe_tax_gen_p0_filterdp.csv"))
write.csv(otu_table(ps_NorthA_gen), file = file.path(path, "NorthA_otu_gen_p0_filterdp.csv"))
write.csv(tax_table(ps_NorthA_gen), file = file.path(path, "NorthA_tax_gen_p0_filterdp.csv"))

write.csv(otu_table(ps_Asia_gen_0.01), file = file.path(path, "Asia_otu_gen_p0.01_filterdp.csv"))
write.csv(tax_table(ps_Asia_gen_0.01), file = file.path(path, "Asia_tax_gen_p0.01_filterdp.csv"))
write.csv(otu_table(ps_Africa_gen_0.01), file = file.path(path, "Africa_otu_gen_p0.01_filterdp.csv"))
write.csv(tax_table(ps_Africa_gen_0.01), file = file.path(path, "Africa_tax_gen_p0.01_filterdp.csv"))
write.csv(otu_table(ps_Europe_gen_0.01), file = file.path(path, "Europe_otu_gen_p0.01_filterdp.csv"))
write.csv(tax_table(ps_Europe_gen_0.01), file = file.path(path, "Europe_tax_gen_p0.01_filterdp.csv"))
write.csv(otu_table(ps_NorthA_gen_0.01), file = file.path(path, "NorthA_otu_gen_p0.01_filterdp.csv"))
write.csv(tax_table(ps_NorthA_gen_0.01), file = file.path(path, "NorthA_tax_gen_p0.01_filterdp.csv"))




# --- Count shared and unique genera across regions ---
gen_Africa <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/Africa_otu_gen_p0.01_filterdp.csv", row.names = 1)
gen_Asia <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/Asia_otu_gen_p0.01_filterdp.csv", row.names = 1)
gen_Europe <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/Europe_otu_gen_p0.01_filterdp.csv", row.names = 1)
gen_NorthA <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/NorthA_otu_gen_p0.01_filterdp.csv", row.names = 1)

gen_Africa <- row.names(gen_Africa)
gen_Asia <- row.names(gen_Asia)
gen_Europe <- row.names(gen_Europe)
gen_NorthA <- row.names(gen_NorthA)

# Create a named list of four genus vectors
genus_list <- list(
  Africa = gen_Africa,
  Asia = gen_Asia,
  Europe = gen_Europe,
  NorthA = gen_NorthA
)

# Compute union
all_genus <- sort(unique(unlist(genus_list)))
all_genus

# Build a presence/absence matrix
presence_matrix <- map_dfc(genus_list, ~ as.integer(all_genus %in% .x)) %>%
  mutate(Genus = all_genus) %>%
  relocate(Genus)

# Tidy table
presence_matrix <- as.data.frame(presence_matrix)
presence_matrix <- presence_matrix %>% column_to_rownames("Genus")
rownames(presence_matrix) <- sub("^[a-z]{1,2}__", "", rownames(presence_matrix))
presence_matrix$Sum <- rowSums(presence_matrix)
presence_matrix <- presence_matrix[order(-presence_matrix$Sum), ]

head(presence_matrix)
table(presence_matrix$Sum)

# Write output
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/compare"
write.csv(presence_matrix, file.path(path, "overlap_gen_0.01_filterdp.csv"))



# --- Differential taxa across continents ---
setwd("data_for_analysis/")
Africa <- read.csv("location_specific/Africa_otu_gen_p0.01_filterdp.csv", row.names = 1)
Asia   <- read.csv("location_specific/Asia_otu_gen_p0.01_filterdp.csv", row.names = 1)
Europe <- read.csv("location_specific/Europe_otu_gen_p0.01_filterdp.csv", row.names = 1)
NorthA <- read.csv("location_specific/NorthA_otu_gen_p0.01_filterdp.csv", row.names = 1)

uniq_prevelance <-
  data.frame(Genus = c(rownames(Africa), rownames(Asia), rownames(Europe), rownames(NorthA)),
             Continent = c(rep("Africa", nrow(Africa)), rep("Asia", nrow(Asia)), rep("Europe", nrow(Europe)), rep("North America", nrow(NorthA)))
             ) %>%
  filter(!Genus %in% c("Other", "Unknown"))
uniq_summary <- table(uniq_prevelance$Genus) %>% data.frame() %>% filter(Freq == 4)

genus_abun <- read.csv("otu_gen_p0.01_filterdp.csv", sep = ",", row.names = 1)
relative_abundance <- read.csv("abundance_prevalence/genus_p0.01_relative_abundance.csv", row.names = 1)
genus_clr <- read.csv("otu_gen_P0.01_filterdp_clr.csv", sep = ",", header = T, row.names = 1)
metadata <- read.csv("ID_meta_final.csv", row.names = 1)

otu <- genus_clr %>%
  rownames_to_column(var = "Taxa") %>%
  pivot_longer(cols = all_of(colnames(genus_clr)), names_to = "ID", values_to = "Abundance") %>%
  left_join(metadata, by = "ID")
continent_list <- c("Africa", "Asia", "Europe", "North America")
taxa_list <- uniq_summary$Var1

# run model (continent pairs)
lmm_result_total <- data.frame()
for(reference in continent_list){
  variable0 = reference
  lmm_result <- data.frame()
  for(taxa in taxa_list){
    df1 <- otu %>% filter(Taxa == taxa)
    df1$Continent <- factor(df1$Continent, levels = c(variable0, setdiff(continent_list, variable0))) # set levels 1 as reference variable
    formula1 <- as.formula("Abundance ~ Continent +  (1 | PRJ) ")
    model <- lmer(formula1, data = df1)
    if (isSingular(model)){check = "Singular"} else { check = "normal"} # check model singular or not
    result_df <- summary(model)$coefficients %>% as.data.frame() %>%
      rownames_to_column(var = "Variable") %>%
      rename_with(~c("Variable", "Estimate", "Std_Error", "df", "t.value",  "p.value")) %>%
      mutate(Lower_CI = Estimate - 1.96 * Std_Error,
             Upper_CI = Estimate + 1.96 * Std_Error,
             Taxa = taxa,
             Reference_variable = variable0,
             Check_model = check
      )
    lmm_result <- rbind(lmm_result, result_df)
  }
  lmm_result_total <- rbind(lmm_result_total, lmm_result)
}
write.table(lmm_result_total, "differential_analysis_res/diff_continent_linear_mix_res.tsv", sep = "\t", row.names = F, quote = F)

res <- lmm_result_total %>%
  mutate(p.value = case_when(Check_model == "Singular" ~ NA, T ~ p.value)) %>%
  filter(!Variable == "(Intercept)" & !Taxa %in% c("Other", "Unknown")) %>%
  group_by(Reference_variable) %>%
  mutate(adjust.p.value = p.adjust(p.value, method = "fdr")) %>%
  arrange(Reference_variable, Variable, Taxa) %>%
  ungroup() %>%
  mutate(Variable = str_replace_all(Variable, "Continent", "")) %>%
  mutate(Significance1 = case_when(adjust.p.value < 0.01 ~ "Sig", T ~ "Non-sig")) %>%
  mutate(Significance2 = case_when(Estimate > 0 & Significance1 == "Sig" ~ "Higher",
                                   Estimate < 0 & Significance1 == "Sig" ~ "Lower",
                                   T ~ "Non-sig"))
write.csv(res, "differential_analysis_res/plotdata_continent_diff_genus.csv")


# --- Abundance and prevalence metrics of genera shared across continents ---
setwd("data_for_analysis")

genus_abun <- read.csv("otu_gen_p0.01_filterdp.csv", sep = ",", row.names = 1)
metadata <- read.csv("ID_meta_final.csv", sep = ",", header = T, row.names = 1)

total_counts <- metadata %>% select(ID, Depth_microfungi) %>% column_to_rownames(var = "ID")
total_counts <- total_counts[colnames(genus_abun), 1]
relative_abundance <- sweep(genus_abun, 2, total_counts, FUN = "/")
write.csv(relative_abundance, "abundance_prevalence/genus_p0.01_relative_abundance.csv")

min(relative_abundance[relative_abundance != 0]) #check min value (except 0)
mean_abun_countinent <- log10(relative_abundance+0.000001) %>% t() %>% data.frame() %>%
  mutate(group = metadata$Continent[match(rownames(.), metadata$ID)]) %>%
  filter(!is.na(group)) %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  data.frame()  %>%
  column_to_rownames(var = "group") %>%
  t()
head(mean_abun_countinent)
write.csv(mean_abun_countinent, "abundance_prevalence/genus_p0.01_log10relative_abundance_mean_countinent.csv")


# prevalence calculate
Africa <- read.csv("location_specific/Africa_otu_gen_p0_filterdp.csv", row.names = 1)
Asia   <- read.csv("location_specific/Asia_otu_gen_p0.01_filterdp.csv", row.names = 1)
Europe <- read.csv("location_specific/Europe_otu_gen_p0.01_filterdp.csv", row.names = 1)
NorthA <- read.csv("location_specific/NorthA_otu_gen_p0.01_filterdp.csv", row.names = 1)


df1 <- list()
df1[[1]] <- Africa
df1[[2]] <- Asia
df1[[3]] <- Europe
df1[[4]] <- NorthA
names(df1) <- c("Africa", "Asia", "Europe", "North America")

taxa_list <- rownames(genus_abun) %>% setdiff(c("Other", "Unknown"))
core_pre <- data.frame()
for(i in 1:4){
  data0 <- df1[[i]] %>% mutate(Pre = rowSums(. > 0)/ncol(.)) %>% rownames_to_column(var = "Genus") %>%
    select(Genus, Pre) %>% filter(Genus %in% taxa_list) %>%
    mutate(Continent = names(df1)[i])
  core_pre <- rbind(core_pre, data0)
}
write.csv(core_pre, "abundance_prevalence/genus_p0.01_prevelance.csv")




# --- Correlation network computation ---
# Executed by 02_run_fastspar.sh



# --- Initial processing of network results ---
location <- "All"
location <- "Europe"
location <- "Africa"
location <- "NorthA"
location <- "Asia"

path_in <- "/data/wangxinyu/ITS_Public/Project/Analysis/fastspar"

cor <- read.delim(file.path(path_in, paste0(location, "_gen_0.01"), paste0(location, "_gen_0.01_cor.tsv")), row.names = 1)
pval <- read.delim(file.path(path_in, paste0(location, "_gen_0.01"), paste0(location, "_gen_0.01_pval.tsv")), row.names = 1)

cor[1:3, 1:3]
pval[1:3, 1:3]

# Build edge list
edge <- cor %>%
  rownames_to_column("taxon1") %>%
  pivot_longer(-taxon1, names_to = "taxon2", values_to = "cor") %>%
  left_join(
    pval %>%
      rownames_to_column("taxon1") %>%
      pivot_longer(-taxon1, names_to = "taxon2", values_to = "pval"),
    by = c("taxon1", "taxon2")
  ) %>%
  filter(taxon1 != taxon2) %>%
  mutate(
    node1 = pmin(taxon1, taxon2), # Normalize taxon order (alphabetical)
    node2 = pmax(taxon1, taxon2)
  ) %>%
  select(node1, node2, cor, pval) %>%
  distinct() # Remove duplicates

# Drop taxonomic prefixes
edge$node1 <- sub("^[a-z]{1,2}__", "", edge$node1)
edge$node2 <- sub("^[a-z]{1,2}__", "", edge$node2)

nrow(edge) # (((num of taxa) * (num of taxa)) - num of taxa)/2
head(edge)

path <- file.path(path_in, paste0(location, "_gen_0.01"))
path
write.csv(edge, file.path(path, paste0(location, "_gen_0.01_edge.csv")))



# --- Shared vs. region-specific significant edges ---
edge_Africa <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/fastspar/Africa_gen_0.01/Africa_gen_0.01_edge.csv", row.names = 1)
edge_Asia <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/fastspar/Asia_gen_0.01/Asia_gen_0.01_edge.csv", row.names = 1)
edge_Europe <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/fastspar/Europe_gen_0.01/Europe_gen_0.01_edge.csv", row.names = 1)
edge_NorthA <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/fastspar/NorthA_gen_0.01/NorthA_gen_0.01_edge.csv", row.names = 1)

# Remove rows with 'Other'/'Unknown', filter by p-value (and |r|), add edge pair key
filter_edges <- function(df, threshold_p = 0.001, threshold_r = 0, prefix = NULL) {
  df %>%
    filter(
      !str_detect(node1, "Other|Unknown") &
      !str_detect(node2, "Other|Unknown") &
      pval <= threshold_p &
      abs(cor) >= threshold_r
    ) %>%
    mutate(
      node1 = pmin(node1, node2),
      node2 = pmax(node1, node2),
      edge = paste0(node1, "-", node2)
    ) %>%
    rename_with(
      ~ if (!is.null(prefix)) paste0(prefix, "_", .x) else .x,
      .cols = "cor"
    )
}

# Filter
edge_Africa_filter <- filter_edges(edge_Africa, prefix = "Africa")
edge_Asia_filter <- filter_edges(edge_Asia, prefix = "Asia")
edge_Europe_filter <- filter_edges(edge_Europe, prefix = "Europe")
edge_NorthA_filter <- filter_edges(edge_NorthA, prefix = "NorthA")

head(edge_Africa_filter)

# Get all non-redundant edges
edge_list <- unique(c(edge_Africa_filter$edge,
                    edge_Asia_filter$edge,
                    edge_Europe_filter$edge,
                    edge_NorthA_filter$edge))
edge_list

# Build a data frame
edge_list <- data.frame(edge = edge_list)
head(edge_list)
nrow(edge_list)

# Merge correlations from four regions
edge_merged <- edge_list %>%
  left_join(select(edge_Africa_filter, edge, Africa_cor), by = "edge") %>%
  left_join(select(edge_Asia_filter, edge, Asia_cor), by = "edge") %>%
  left_join(select(edge_Europe_filter, edge, Europe_cor), by = "edge") %>%
  left_join(select(edge_NorthA_filter, edge, NorthA_cor), by = "edge")

head(edge_merged)
nrow(edge_merged)

# Count positive/negative counts
edge_merged2 <- edge_merged %>%
  rowwise() %>%
  mutate(
    pos_count = sum(c_across(ends_with("_cor")) > 0, na.rm = TRUE),
    neg_count = sum(c_across(ends_with("_cor")) < 0, na.rm = TRUE),
    all_count = pos_count + neg_count,
    direction = case_when(
      pos_count > 0 & neg_count == 0 ~ "positive",
      neg_count > 0 & pos_count == 0 ~ "negative",
      pos_count > 0 & neg_count > 0  ~ "mixed",
      TRUE ~ "undetermined"
    )
  ) %>%
  ungroup()

# Assess sharing per edge across regions
# For each edge, set shared = 1 if any other region has the same correlation sign
edge_merged3 <- edge_merged2 %>%
  rowwise() %>%
  mutate(
    # Presence (non-NA)
    Africa_presence = ifelse(!is.na(Africa_cor), 1, 0),
    Asia_presence   = ifelse(!is.na(Asia_cor),   1, 0),
    Europe_presence = ifelse(!is.na(Europe_cor), 1, 0),
    NorthA_presence = ifelse(!is.na(NorthA_cor), 1, 0),

    # Presence with consistent sign to any other region
    Africa_shared = ifelse(
      !is.na(Africa_cor) &&
      any(sign(Africa_cor) == sign(c(Asia_cor, Europe_cor, NorthA_cor)), na.rm = TRUE) &&
      sign(Africa_cor) != 0,
      1, 0
    ),
    Asia_shared = ifelse(
      !is.na(Asia_cor) &&
      any(sign(Asia_cor) == sign(c(Africa_cor, Europe_cor, NorthA_cor)), na.rm = TRUE) &&
      sign(Asia_cor) != 0,
      1, 0
    ),
    Europe_shared = ifelse(
      !is.na(Europe_cor) &&
      any(sign(Europe_cor) == sign(c(Africa_cor, Asia_cor, NorthA_cor)), na.rm = TRUE) &&
      sign(Europe_cor) != 0,
      1, 0
    ),
    NorthA_shared = ifelse(
      !is.na(NorthA_cor) &&
      any(sign(NorthA_cor) == sign(c(Africa_cor, Asia_cor, Europe_cor)), na.rm = TRUE) &&
      sign(NorthA_cor) != 0,
      1, 0
    )
  ) %>%
  ungroup()

edge_merged3 <- as.data.frame(edge_merged3)
edge_merged3 <- edge_merged3[order(-edge_merged3$all_count), ]

head(edge_merged3)

table(edge_merged3$all_count)
table(edge_merged3$direction)
table(edge_merged3$pos_count)
table(edge_merged3$neg_count)



path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/compare"
write.csv(edge_merged3, file.path(path, "overlap_edge_0.01_filterdp.csv"))
