
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
                                 simple_min_max = F) {
  if (dir.exists(root)) warning('Root folder already exists!')

  dir.create(root)
  dir.create(paste0(root, '/new'))

  FUN <- function(x) {
    a <- cnvs[x]

    pt <- paste0(root, '/new/samp', a$sample_ID, '_st', a$start, '.png')

    if (!file.exists(pt)) {

      dt <- plot_cnv(a, samps[sample_ID == a[, sample_ID], ], snps = snps,
                     shrink_lrr = shrink_lrr, simple_min_max = simple_min_max)

      if (nrow(dt) == 0) {
        warning('no image saved for cnv: ', a)
        return(data.table())
      }

      imager::save.image(imager::as.cimg(dt), pt)
    }

    # to be tested, might be unstable when called by multiple workers
    if (x %% 100 == 0) gc()
  }

  # save images using BiocParallel
  null <- BiocParallel::bplapply(1:nrow(cnvs), FUN)

}
