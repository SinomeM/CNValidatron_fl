library(data.table)
library(ggplot2)

setwd("~/Downloads/CNValidatron_fl")
cnv <- fread("./data/put_cnv.txt")
tmpd <- QCtreeCNV:::get_region_tabix(cnv$chr, cnv$start, cnv$end,
                                    "./data/toy_snp_array_noCNV.txt.gz")
tmpc <- fread("./data/put_cnv.txt")
tmpl <- fread("./data/loci.txt")

source("./R/create_image.R")
tmp <- cnv_image_matrix(dt = tmpd, region = tmpl, n = 32, eps = 2,
                        sides = 4, adj = T, path_img = "./tmp/prova.png")
tmp
