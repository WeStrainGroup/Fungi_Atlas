# Purpose: Summarize benchmark accuracy metrics as means and standard deviations for reporting.


library(tidyverse)

taxa <- read.csv("/home/wangxy/ITS_benchmark/analysis/taxa_stats.csv", row.names = 1)
head(taxa)

taxa %>%
  filter(method == "mycogap") %>%
  summarize(
    across(
      c(Precision, Recall, F1_score),
      list(
        mean = ~round(mean(.x, na.rm = TRUE), 3),
        sd   = ~round(sd(.x, na.rm = TRUE), 3)
      ),
      .names = "{.col}_{.fn}"
    )
  )

com <- read.csv("/home/wangxy/ITS_benchmark/analysis/com_stats.csv", row.names = 1)
head(com)


# Calculate mean and sd for mycogap-suffixed columns
com %>%
  summarize(
    across(
      contains("mycogap"),
      list(
        mean = ~round(mean(.x, na.rm = TRUE), 3),
        sd   = ~round(sd(.x, na.rm = TRUE), 3)
      ),
      .names = "{.col}_{.fn}"
    )
  )
