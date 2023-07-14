# install latest QCtreeCNV
if (F) devtools::install_github("sinomem/QCtreeCNV")

setwd('~/Documents/CNValidatron_fl')

# libs
library(data.table)

# Load CNV and samples table from UKB export
cnvs <- fread('~/Documents/UKB_data/edited/cvns.txt')
samples <- fread('~/Documents/UKB_data/edited/samples.txt')
loci <- fread('~/Documents/UKB_data/edited/loci.txt')

cnvs <- cnvs[numsnp >= 20,]

# check everything is fine
all(cnvs$sample_ID %in% samples$sample_ID)

# R torch
library(torch)
cuda_is_available()

# test the function to create the PNG plots of the CNV calls
source('./R/create_image.R')
a <- cnvs[sample(1:nrow(cnvs), 1)]
a
int_data <- QCtreeCNV:::get_region_tabix(a$chr, 0, 1000000000,
                                         samples[sample_ID == a$sample_ID, file_path_tabix] )

cnv_image_matrix(dt = int_data, region = a, n = 50, eps = 4, sides = 0.5,
                 adj = T, path_img = 'tmp/prova3.png')
