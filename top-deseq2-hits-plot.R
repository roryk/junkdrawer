library(tidyr)
library(ggplot2)
library(dplyr)
ncounts = counts(dds, normalized=TRUE)
top20 = head(subset(p1_treatment_df, mgi_symbol != "" & !is.na(mgi_symbol)), 20)
top20 = top20[, c("mgi_symbol", "ensembl_gene_id")]
top20 = top20 %>%
  left_join(as.data.frame(ncounts) %>%
            tibble::rownames_to_column("ensembl_gene_id")) %>%
  gather(Name, count, -mgi_symbol, -ensembl_gene_id) %>%
  left_join(summarydata)

ggplot(top20, aes(Name, count, fill=agetreat)) +
  geom_bar(stat='identity') +
  facet_wrap(~mgi_symbol, scale="free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
