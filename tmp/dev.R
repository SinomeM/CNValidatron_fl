
devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

devtools::install()


library(data.table)

cnvs <- fread('../../UKB_GW_CNVs/cnvs_pred.txt')[pred %in% 2:3 & pred_prob >= 0.9, ]
# mild filters on outliers
cnvs <- cnvs[length <= 10000000, ]
cnvs <- cnvs[!sample_ID %in% cnvs[ , .N, by = sample_ID][N >10, sample_ID], ]
fwrite(cnvs[, .(chr, start, end, GT)], '../../UKB_GW_CNVs/igv/cnvs.bed', col.names = F, sep = '\t')

devtools::load_all()
dt <- cnvrs_iou(cnvs, QCtreeCNV::hg19_chr_arms)
dt
cnvs <- dt[[1]]
cnvrs <- dt[[2]]

