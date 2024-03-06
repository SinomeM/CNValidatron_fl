
devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

devtools::install()

library(data.table)
pt <- '~/Documents/CNValidatron_fl/tmp/prova'
samps <- fread('~/Documents/UKB_data/edited/samples.txt')[sample(200),]
cnvs <- fread('~/Documents/CNValidatron_trained_models/all_visual_inspection_min40snps.txt')[
          sample_ID %in% samps$sample_ID,]
snps <- fread('~/Documents/UKB_data/snppos_filtered.txt')

save_pngs_prediction(pt, cnvs, samps, snps,  w = 64)


devtools::load_all()
pred <- make_predictions(luz::luz_load('~/Documents/CNValidatron_trained_models/fitted_dropout_all.rds'),
                         pt, cnvs)
pred


library(data.table)
samps <- fread('~/Documents/CNValidatron_trained_models/second_iteration/all_samples_15k.txt')
cnvs <- fread('~/Documents/CNValidatron_trained_models/second_iteration/all_training_example_8k.txt')
a <- cnvs[1]
b <- samps[sample_ID == a$sample_ID, ]
plot_cnv(a, b)



library(data.table)
cnvs <- fread('../../UKB_GW_CNVs/cnvs_with_cnvrs.txt')

devtools::load_all()
# debugonce(cnvrs_iou)
dt <- cnvrs_iou(cnvs, QCtreeCNV::hg19_chr_arms)
dt

devtools::load_all()
a <- get_igraph_objs(dt, 0.75, 1, 99, 'arm1')
a

dt
min_iou <- 0.75
leiden_res <- 1
ii <- 11
arm <- 'arm1'

get_igraph_objs <- function(dt, min_iou, leiden_res, ii, arm, type = '') {

  # get all overlapping CNVs
  setkey(dt, start, end)
  dt[, cix := 1:.N]
  dt_s <- foverlaps(dt[, .(cix, start, end)], dt[, .(cix, start, end)])
  colnames(dt_s) <- c('cnvA', 'stA', 'enA', 'cnvB', 'stB', 'enB')
  dt_s <- dt_s[cnvA != cnvB, ]

  dt_s[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
                (pmax(enA, enB) - pmin(stA, stB))]

  # cnvs with connections
  dt_g <- dt_s[iou >= min_iou, ]
  setorder(dt_g, cnvA)
  dt_r <- dt[cix %in% dt_g[, unique(c(cnvA, cnvB))], ]
  # cnvs without connections
  dt_m <- dt[!cix %in% dt_g[, unique(c(cnvA, cnvB))], ]

  # create the igraph network
  g <- graph_from_data_frame(dt_g[, .(cnvA, cnvB)], directed = F)
  gr <- cluster_leiden(g, resolution_parameter = leiden_res)

  # merge dt_r and igraph object here, careful, not all CNVs that got in are necessarily in dt_r
  dt_r[, CNVR := paste0(arm, '_', ii, '_', ifelse(type != '', paste0(type, '_'), ''),
                        membership(gr))]
  if (dt_m[, .N] > 0) {
    dt_s[, CNVR := paste0(arm, '_', ii, '_singleton_', 1:.N)]
    dt <- rbind(dt_s, dt_r)
    return(list(g, gr, dt))
  }

  return(list(g, gr, dt_r))
}
