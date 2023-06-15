library(data.table)
library(dplyr)
library(GenomicRanges)

gr <- rtracklayer::import("frank_perturb.bed")
bam <- rtracklayer::import("GEX_cd69p1_rna.bam")

data.frame(
  gr,
  reads = countOverlaps(gr, bam)
) %>% mutate(igv_lookup = paste0(seqnames, ":", start, "-", end)) %>%
  arrange(reads)
