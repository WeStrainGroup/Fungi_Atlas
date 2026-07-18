# Purpose: Plot the genera most strongly correlated with the MDS2 ordination axis.
# Panels: Figure S7E.

library(tidyverse)
library(ggplot2)

# --- Plot genera associated with the MDS2 axis ---
res_cor_mds2_taxa <- read.csv("sensitive_analysis/mds2_driver_genus.csv")
p_driver_feature2 <- ggplot(res_cor_mds2_taxa, aes(y = Genus_Unique, x = r_value, fill = cor_direction)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = round(r_value, 2),
                hjust = ifelse(r_value > 0, -0.1, 1.1)),
            color = "black", size = 3) +
  theme_bw() +
  scale_fill_manual(values = c("neg" = "#6a91b8", "pos" = "#b3656e")) +
  xlim(-1, 1) +
  labs(x = "Spearman's r", y = "") +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(face = "italic")
  ) +
  scale_y_discrete(labels = function(x) gsub(".*__(.*)", "\\1", x))
p_driver_feature2
ggsave(plot = p_driver_feature2, "plot/mds2_driver_genera_top10.pdf", width = 3, height = 4)
