## (RE)START FROM SCRATCH THE DESIGN ##

setwd('~/Documents/CNValidatron_fl')
library(data.table)

# Load CNV and samples table from UKB export
cnvs <- fread('~/Documents/UKB_data/edited/cvns.txt')
samples <- fread('~/Documents/UKB_data/edited/samples.txt')
loci <- fread('~/Documents/UKB_data/edited/loci.txt')
snps <- fread('~/Documents/UKB_data/snppos_filtered.txt')
cnvs <- cnvs[numsnp >= 20,]

# load visual validations
tmp <- list.files('~/Documents/UKB_data/visual_inspection2/')
tmp <- tmp[grep('visual', tmp)]
vi_cnv <- data.table()
for (i in tmp)
  vi_cnv <- rbind(vi_cnv,
                  fread(paste0('~/Documents/UKB_data/visual_inspection2/', i)))

# check for duplicates
vi_cnv[, ':=' (V1 = NULL, index = NULL)]
nrow(unique(vi_cnv)) == nrow(vi_cnv)
nrow(unique(vi_cnv[, .(sample_ID, chr, start, GT)])) == nrow(vi_cnv)
tmp <- vi_cnv[duplicated(vi_cnv[, .(sample_ID, chr, start, GT)])]
vi_cnv[sample_ID %in% tmp$sample_ID & chr %in% tmp$chr & start %in% tmp$start &
       GT %in% tmp$GT, ]

# remember to deal with the modified boundaries
vi_cnv[Visual_Output == -8, ]

# It is useful that here each locus is unique
length(unique(cnvs$locus)) == nrow(cnvs)


## STEP 1 ##

source('./R/02_new_dev_functions.R')

tmp <- load_snps_tbx(cnvs[1], samples[sample_ID == cnvs[1, sample_ID], ], snps = snps)


tmp <- plot_cnv(cnvs[1], samples[sample_ID == cnvs[1, sample_ID], ], snps = snps)
tmp
