
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
min_iou <- 0.75
leiden_res <- 1
amr <- 'adas'
ii <- 11

merge_cnvrs <- function(cnvrs, cnvs, min_iou, leiden_res, arm, ii) {

  # compute iou between all pairs of overlapping CNVRs
  setkey(cnvrs, start, end)
  dt <- foverlaps(cnvrs[, .(start, end, CNVR)], cnvrs[, .(start, end, CNVR)])
  setnames(dt, c('CNVR', 'start', 'end', 'i.CNVR', 'i.start', 'i.end'),
           c('cnvrA', 'stA', 'enA', 'cnvrB', 'stB', 'enB'))
  dt <- dt[cnvrA != cnvrB & !(is.na(stA) | is.na(stB)), ]
  dt[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
              (pmax(enA, enB) - pmin(stA, stB))]

  if (dt[, .N] == 0) return(list(cnvs, cnvrs))

  message('Merging come CNVRs!')

  # I can treat CNVRs as CNVs and use the same idea of the get_igraph_objs() function

  # cnvrs with connections
  dt_g <- dt[iou >= min_iou, ]
  setorder(dt_g, cnvrA)
  dt_r <- cnvrs[CNVR %in% dt_g[, unique(c(cnvrA, cnvrB))], ]
  dt_r[, CNVR_old := CNVR]
  # cnvs without connections
  dt_m <- cnvrs[!CNVR %in% dt_g[, unique(c(cnvrA, cnvrB))], ]

  # create the igraph network
  g <- graph_from_data_frame(dt_g[, .(cnvrA, cnvrB)], directed = F)
  gr <- cluster_leiden(g, resolution_parameter = leiden_res)
  a <- membership(gr)
  a <- data.table(CNVR = names(a), gr = a)
  # re construct the CNVRs table
  dt_r <- merge(dt_r, a, by = 'CNVR')
  dt_r[, CNVR := paste0(arm, '_', ii, '_merge_', gr)]
  dt_r[, gr := NULL]

  # Updated CNVRs info in the CNVs table
  for (i in dt_r[, unique(CNVR)])
    cnvs[CNVR %in% dt_r[CNVR == i, CNVR_old], CNVR := i]

  # recreate CNVRs table is from updated CNVs table
  cnvrs <- create_cnvrs(cnvs)

  return(list(cnvs, cnvrs))
}
