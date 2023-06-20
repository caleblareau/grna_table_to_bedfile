library(data.table)
library(dplyr)
library(GenomicRanges)
library(GenomicAlignments)
library(stringr)
library(BuenColors)
# Import data and set up
bed_gr <- rtracklayer::import("frank_perturb.bed")
frank_anno <- read.table("cd69p_barcodes-by-sgrna-call.csv", sep = ",", header = TRUE) %>%
  mutate(barcode = gsub("cd69p1_", "", X)) %>%
  mutate(target = str_split_fixed(sgrna.class, "-", 2)[,1])
bam <-GenomicAlignments::readGAlignments("GEX_cd69p1_rna.bam", param = ScanBamParam(tag=c("CB")))

# Now do a bunch of loops
lapply(1:length(bed_gr), function(i){
  bam_ss <- bam[findOverlaps(bam, bed_gr[i])@from]
  name <- bed_gr$name[i]
  target_barcodes <- frank_anno %>% filter(target == name) %>%
    pull(barcode)
  non_target_barcodes <- frank_anno %>% filter(target != name) %>%
    pull(barcode)
  rle_target <- coverage(bam_ss[bam_ss@elementMetadata$CB %in% target_barcodes])[[as.character(bed_gr@seqnames[i])]]
  rle_nontarget <- coverage(bam_ss[bam_ss@elementMetadata$CB %in% non_target_barcodes])[[as.character(bed_gr@seqnames[i])]]
  
  # Subset coordinates
  range <- bed_gr@ranges[i]@start:(bed_gr@ranges[i]@start + bed_gr@ranges[i]@width)
  as.numeric(rle_target[range])
  as.numeric(rle_nontarget[range])
  data.frame(range, 
             target =   as.numeric(rle_target[range]),
             offtarget =  as.numeric(rle_nontarget[range])
             ) %>%
    reshape2::melt(id.var = "range") %>%
    ggplot(aes(x = range, y = value, color = variable)) + 
    geom_line() + pretty_plot() +
    facet_wrap(~variable, scales = "free_y") + ggtitle(name, paste0(": ", paste0(bed_gr[i], collapse = "_"))) +
    labs(x = "Genomic Coordinate", y = "coverage per base pair") -> p1
  cowplot::ggsave2(p1, file = paste0("plots/",name, "_", as.character(i), ".png"), width = 6, height = 2.5)
})