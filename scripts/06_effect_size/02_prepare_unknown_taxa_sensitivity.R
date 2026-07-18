# Purpose: Remove the aggregated Unknown genus and rebuild CLR profiles and robust Aitchison distances for sensitivity tests.
# Panels: Figure S7B; also prepares the drop-Unknown CLR input required for Figure S7G.

library(data.table)
library(vegan)



otu_gen <- fread("/data/wangxinyu/Fungi_Atlas/Analysis/Main/data/genus/otu_gen_p0.001_filterdp.csv", data.table = F, header = TRUE)
otu_gen <- otu_gen %>% column_to_rownames("V1") %>% t() %>% as.data.frame()
otu_gen[1:5, 1:5]
dim(otu_gen) #  37417   688

otu_gen_dropunknown <- otu_gen %>% select(-Unknown)
dim(otu_gen_dropunknown) # 37417   687


dist_gen_dropunknown <- vegdist(otu_gen_dropunknown, MARGIN = 1, method = "robust.aitchison")
str(dist_gen_dropunknown)


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/unknown_check"
saveRDS(dist_gen_dropunknown, file.path(path, "dist_gen_dropunknown.rds"))





otu_gen <- fread("/data/wangxinyu/Fungi_Atlas/Analysis/Main/data/genus/otu_gen_p0.01_filterdp.csv", data.table = F, header = TRUE)
otu_gen <- otu_gen %>% column_to_rownames("V1") %>% t() %>% as.data.frame()
otu_gen[1:5, 1:5]
dim(otu_gen) #  37417  200


otu_gen_dropunknown <- otu_gen %>% select(-Unknown)
dim(otu_gen_dropunknown) # 37417  199
otu_gen_dropunknown[1:5, 1:5]


otu_gen_dropunknown_clr <- vegan::decostand(otu_gen_dropunknown, MARGIN = 1, method = "clr", pseudocount = 1)
otu_gen_dropunknown_clr[1:5, 1:5]


otu_gen_dropunknown_clr[1:5, 1:5]
rowSums(otu_gen_dropunknown_clr[1:10, ])
apply(otu_gen_dropunknown_clr[1:10, ], 1, sd)


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/unknown_check"
write.csv(otu_gen_dropunknown_clr, file.path(path, "otu_gen_dropunknown_clr.csv"))
