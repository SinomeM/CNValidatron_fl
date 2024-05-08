
#' Save the PNG Plots of CNVs for Prediction
#'
#' This function is used to create the dataset of PNG images for prediction
#'
#' @param root root folder for the dataset. Must not exists.
#' @param cnvs cnv data.table in the usual format
#' @param samps sample list in usual format
#' @param snps snps in the usual format
#' @param shrink_lrr see load_snps_tbx(), should be the same the model used in training
#'
#' @export
#'
#' @import data.table

save_pngs_prediction <- function(root, cnvs, samps, snps, shrink_lrr = 0.2,
                                 simple_min_max = F, no_parall = F, batches = 1000) {
  if (dir.exists(root)) warning('Root folder already exists!')

  dir.create(root)
  #dir.create(paste0(root, '/new'))

  if (!'batch' %in% colnames(cnvs)) cnvs[, batch := sample(1:batches, .N, replace = T)]
  for (i in 1:batches) dir.create(paste0(root, '/batch', i), showWarnings = F)

  FUN <- function(x) {
    a <- cnvs[x]

    dt <- plot_cnv(a, samps[sample_ID == a[, sample_ID], ], snps = snps,
                   shrink_lrr = shrink_lrr, simple_min_max = simple_min_max)
    n_real_snps <- dt[[2]]
    dt <- dt[[1]]

    pt <- paste0(root, '/new/batch', a$batch, '/samp', a$sample_ID, '_st', a$start,
                 '_nsnp', n_real_snps, '.png')

    #if (!file.exists(pt)) {
    #  if (nrow(dt) == 0) {
    #    warning('no image saved for cnv: ', a)
    #    return(data.table())
    #}

    imager::save.image(imager::as.cimg(dt), pt)
    #}

    if (x %% 100 == 0) gc()
  }

  if (no_parall)  null <- lapply(1:nrow(cnvs), FUN)
  else null <- BiocParallel::bplapply(1:nrow(cnvs), FUN)
  # save images using BiocParallel

}
