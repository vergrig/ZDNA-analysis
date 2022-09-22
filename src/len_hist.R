library('ggplot2')

DATA_DIR <- "./"
OUT_DIR <- "./results"

#NAME <- 'H3K9me9_SJSA1.ENCFF921OTR.hg38'
#NAME <- 'H3K9me9_SJSA1.ENCFF921OTR.hg19'
#NAME <- 'H3K9me9_SJSA1.ENCFF157SWY.hg38'
#NAME <- 'H3K9me9_SJSA1.ENCFF157SWY.hg19'
#NAME <- 'H3K9me9_SJSA1.intersect_with_DeepZ'
NAME <- 'DeepZ'


bed_df <- read.delim(paste0(DATA_DIR, NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end')
bed_df$len <- bed_df$end - bed_df$start

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.pdf'), path = OUT_DIR)
