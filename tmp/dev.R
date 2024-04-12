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

if (F) {
  library(data.table)
  tdel <- fread("~/Documents/CNValidatron_trained_models/third_iteration/visual_inspection/true_dels_vi_res.txt")
  tdup <- fread("~/Documents/CNValidatron_trained_models/third_iteration/visual_inspection/true_dups_vi_res.txt")
  fals <- fread("~/Documents/CNValidatron_trained_models/third_iteration/visual_inspection/false_vi_res.txt")
  samp <- fread('~/Documents/CNValidatron_trained_models/third_iteration/samples.txt')
  snps <- fread('~/Documents/UKB_data/data/snppos_filtered.txt')

  # small
  a <- tdel[length <=  100000 & vo == 1, ][sample(1:.N, 1), ]
  a <- tdup[length <=  100000 & vo == 1, ][sample(1:.N, 1), ]
  a <- fals[length <=  100000 & vo == 2, ][sample(1:.N, 1), ]
  # medium/large
  a <- tdel[length >  100000 & vo == 1, ][sample(1:.N, 1), ]
  a <- tdup[length >  100000 & vo == 1, ][sample(1:.N, 1), ]
  a <- fals[length >  100000 & vo == 2, ][sample(1:.N, 1), ]
  # large
  a <- tdel[length >  1000000 & vo == 1, ][sample(1:.N, 1), ]
  a <- tdup[length >  1000000 & vo == 1, ][sample(1:.N, 1), ]
  a <- fals[length >  1000000 & vo == 2, ][sample(1:.N, 1), ]

  devtools::load_all()
  plot_cnv(a, samp[sample_ID == a$sample_ID, ], snps, tmp_plot = 3, shrink_lrr = 0.2)
}
