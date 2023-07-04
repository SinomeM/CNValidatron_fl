# install latest QCtreeCNV
if (F) devtools::install_github("sinomem/QCtreeCNV")

# libs
library(data.table)

# Load CNV and samples table from UKB export
cnvs <- fread('~/Documents/UKB_data/export/cnvs.txt')
samples <- fread('~/Documents/UKB_data/export/samples.txt')

# check everything is fine
all(cnvs$sample_ID %in% samples$sample_ID)

# R torch
library(torch)
cuda_is_available()

# test the function to create the PNG plots of the CNV calls
cnvs[, local_fp := paste0('~/Documents/UKB_data/export/tabix/', sample_ID, '.tabix.gz')]
cnvs[, chr := gsub('chr', '', chr)]
a <- cnvs[sample(1:nrow(cnvs), 1)]
int_data <- QCtreeCNV:::get_region_tabix(a$chr, a$start, a$stop, a$local_fp)
# get_region_tabix seems to work

source('./R/create_image.R')
cnv_image_matrix()