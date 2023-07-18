## (RE)START FROM SCRATCH THE DESIGN ##

setwd('~/Documents/CNValidatron_fl')
library(data.table)
library(ggplot2)

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
# produce the plots in a format the algorithm can undrstand best

source('./R/02_new_dev_functions.R')

# option 1: pixelated image of classic BAF LRR in a square figure
ii <- 12
cc <- vi_cnv[numsnp > 40 & Visual_Output == 1, ][ii]

ior <- 3; shlrr <- 0.25; width <- 50

print(check_cnv(cc, samples[sample_ID == cc[, sample_ID], ], snps = snps,
                in_out_ratio = ior, shrink_lrr = shlrr))

tmp <- plot_cnv(cc, samples[sample_ID == cc[, sample_ID], ], snps = snps, tmp_plot = 2,
                w = width, in_out_ratio = ior, shrink_lrr = shlrr)

## STEP 2 ##
# create and store the data for each CNV

# - In torch the plots will be tensors
# - In torch the data handling is done via a dataloader()
# - Since the function torchvision::image_folder_dataset() exists, the
#   easiest thing to do is save all plots as PNG in the expected file
#   structure

# just to test data structure etc really, way too small to even early testing
# of an actual model
first_test_cnvs <- vi_cnv[numsnp > 40]

save_pngs_dataset('/home/simone/Documents/CNValidatron_fl/tmp/first_test',
                  first_test_cnvs, samples, snps)
