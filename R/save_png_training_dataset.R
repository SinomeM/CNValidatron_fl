#' Create Dataset of Validated CNVs
#'
#' This function is used to create the dataset of PNG images for training.
#' At the moment not all options for plot_cnv() are available here.
#'
#' @param root root folder for the dataset. Must not exists.
#' @param cnvs cnv data.table in the usual format plus the column 'vo'
#' @param samps sample list in usual format
#' @param snps snps in the usual format
#' @param shrink_lrr see load_snps_tbx()
#' @param flip_chance probability of saving a flipped example as well
#' @param noise_chance probability of saving a "noised" example as well
#' @param noise_lvl level of noise introduced in the PNG
#'
#' @export
#'
#' @import data.table

save_pngs_dataset <- function(root, cnvs, samps, snps, shrink_lrr = 0.2, flip_chance = 0.5,
                              noise_chance = 0, noise_lvl = 0.1,
                              simple_min_max = F) {
  if (dir.exists(root)) warning('Root folder already exists!')

  dir.create(root)
  dir.create(paste0(root, '/true_del')); dir.create(paste0(root, '/true_dup'))
  dir.create(paste0(root, '/false'))

  FUN <- function(x) {
    a <- cnvs[x]

      dt <- plot_cnv(a, samps[sample_ID == a[, sample_ID], ], snps = snps,
                     shrink_lrr = shrink_lrr, simple_min_max = simple_min_max)
      n_real_snps <- dt[[2]]
      dt <- dt[[1]]

    if (a$GT == 1) {
      if (a$vo == 1) pt <- paste0(root, '/true_del/samp', a$sample_ID,
                                  '_st', a$start, '_nsnp', n_real_snps, '.png')
      if (a$vo == 2) pt <- paste0(root, '/false/samp', a$sample_ID,
                                  '_st', a$start, '_nsnp', n_real_snps, '.png')
    }
    if (a$GT == 2) {
      if (a$vo == 1) pt <- paste0(root, '/true_dup/samp', a$sample_ID,
                                  '_st', a$start, '_nsnp', n_real_snps, '.png')
      if (a$vo == 2) pt <- paste0(root, '/false/samp', a$sample_ID,
                                  '_st', a$start, '_nsnp', n_real_snps, '.png')
    }

    if (!file.exists(pt)) {
      if (nrow(dt) == 0) {
        warning('no image saved for cnv: ', a)
        return(data.table())
      }

      imager::save.image(imager::as.cimg(dt), pt)

      # Data augmentation 1, image flipping
      if (runif(1) <= flip_chance) {
        dt[, x := abs(x-(max(x)+1))] # flip the x axis
        pt <- gsub('\\.png', '_flip\\.png', pt)

      # Data augmentation 2, add some random noise
      if (runif(1) <= noise_chance) {
        nv <- rnorm(w*w, sd = 0.3) * noise_lvl
        dt[, value := value + nv]
        pt <- gsub('\\.png', '_noised\\.png', pt)
      }

        imager::save.image(imager::as.cimg(dt), pt)
      }

      # to be tested, might be unstable when called by multiple workers
      if (x %% 100 == 0) gc()
    }
  }

  # save images using BiocParallel
  null <- BiocParallel::bplapply(1:nrow(cnvs), FUN)

}


