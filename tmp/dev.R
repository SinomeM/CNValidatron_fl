devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

devtools::install()


# test CNVRs
if (F) {
  library(data.table)
  devtools::load_all()
  cnvs <- fread('../../UKB_GW_CNVs/cnvs_pred.txt')[pred %in% 2:3 & pred_prob >= 0.9, ]
  # mild filters on outliers
  cnvs <- cnvs[length <= 10000000, ]
  cnvs <- cnvs[!sample_ID %in% cnvs[ , .N, by = sample_ID][N >10, sample_ID], ]
  fwrite(cnvs[, .(chr, start, end, GT)], '../../UKB_GW_CNVs/igv/cnvs.bed', col.names = F, sep = '\t')

  dt <- cnvrs_iou(cnvs, QCtreeCNV::hg19_chr_arms)
  dt
  cnvs <- dt[[1]]
  cnvrs <- dt[[2]]
}


# test new PNGs
if (F) {
  library(data.table)
  dels <- fread('../../UKB_GW_CNVs/calibration/dels1.txt')
  dups <- fread('../../UKB_GW_CNVs/calibration/dups1.txt')
  samples <- fread('../../UKB_GW_CNVs/calibration/samples.txt')
  # get a random deletion from the calbration set
  a <- dels[numsnp > 50, ][sample(1:.N, 1)]
  b <- samples[sample_ID == a$sample_ID, ]
  devtools::load_all()
  # debugonce(plot_cnv)
  plot_cnv(a, b, tmp_plot = 1, shrink_lrr = 0.2)
  plot_cnv(a, b, tmp_plot = 2, shrink_lrr = 0.2)
  plot_cnv(a, b, tmp_plot = 3, shrink_lrr = 0.2)
  # get a random duplication from the calbration set
  a <- dups[numsnp > 50, ][sample(1:.N, 1)]
  b <- samples[sample_ID == a$sample_ID, ]
  devtools::load_all()
  # debugonce(plot_cnv)
  plot_cnv(a, b, tmp_plot = 1, shrink_lrr = 0.2)
  plot_cnv(a, b, tmp_plot = 2, shrink_lrr = 0.2)
  plot_cnv(a, b, tmp_plot = 3, shrink_lrr = 0.2)
}

