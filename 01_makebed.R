library(data.table)
library(annotables)
library(dplyr)

chr_vec <- grcm38$chr; names(chr_vec) <- grcm38$symbol
trm <- fread("TRM_Perturb_brie_subset_detailed.csv")
trm$chr_basic <- as.numeric(chr_vec[as.character(trm$Target.Gene.Symbol)])

data.frame(
  chr = trm$chr_basic,
  start = trm$Position.of.Base.After.Cut..1.based. -100,
  end = trm$Position.of.Base.After.Cut..1.based. + 100,
  gene = trm$Target.Gene.Symbol
) %>% arrange(chr, start) %>%
  write.table(file = "frank_perturb.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
