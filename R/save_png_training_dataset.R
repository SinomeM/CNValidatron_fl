#' Create Dataset of Validated CNVs
#'
#' This function is used to create the dataset of PNG images for training.
#' At the moment not all options for plot_cnv() are available here.
#'
#' @param root root folder for the dataset. Must not exists.
#' @param cnvs cnv data.table in the usual format
#' @param samps sample list in usual format
#' @param snps snps in the usual format
#' @param w see plot_cnv()
#' @param in_out_ratio see load_snps_tbx()
#' @param shrink_lrr see load_snps_tbx()
#' @param flip_chance probability of saving a flipped example as well
#'
#' @export
#'
#' @import data.table

save_pngs_dataset <- function(root, cnvs, samps, snps, w = 64, in_out_ratio = 3,
                              shrink_lrr = 0.2, flip_chance = 0.5) {
  if (dir.exists(root)) stop('Root folder already exists. Delete existing folder or provide a different path')

  dir.create(root)
  if(5 %in% cnvs$Visual_Output) dir.create(paste0(root, '/not_eval'))
  else {
    dir.create(paste0(root, '/true_del')); dir.create(paste0(root, '/true_dup'))
    dir.create(paste0(root, '/unk_dup')); dir.create(paste0(root, '/unk_del'))
    dir.create(paste0(root, '/false'))
  }


  FUN <- function(x) {
    a <- cnvs[x]

    if (a$Visual_Output == 5) pt <- paste0(root, '/not_eval/', a$sample_ID, '_', a$start, '.png')

    if (a$GT == 1) {
      if (a$Visual_Output == 1) pt <- paste0(root, '/true_del/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 2) pt <- paste0(root, '/false/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 3) pt <- paste0(root, '/unk_del/', a$sample_ID, '_', a$start, '.png')
    }
    if (a$GT == 2) {
      if (a$Visual_Output == 1) pt <- paste0(root, '/true_dup/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 2) pt <- paste0(root, '/false/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 3) pt <- paste0(root, '/unk_dup/', a$sample_ID, '_', a$start, '.png')
    }

    dt <- plot_cnv(a, samps[sample_ID == a[, sample_ID], ], snps = snps,
                   w = w, in_out_ratio = in_out_ratio, shrink_lrr = shrink_lrr)

    if (nrow(dt) == 0) {
      warning('no image saved for cnv: ', a)
      return(data.table())
    }

    dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis
    imager::save.image(imager::as.cimg(dt), pt)

    # Data augmentation 1, image flipping
    if (runif(1) <= flip_chance) {
      dt[, x := abs(x-(max(x)+1))] # flip the x axis
      pt <- gsub('\\.png', '_flip\\.png', pt)
      imager::save.image(imager::as.cimg(dt), pt)
    }
  }

  # save images using BiocParallel
  null <- BiocParallel::bplapply(1:nrow(cnvs), FUN)

}


