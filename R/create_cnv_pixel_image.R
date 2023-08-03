#' Create the matrix for the PNG image
#'
#' This function create the pixel matrix that can be saved
#' as a PNG for further use. All processing is done in the
#' function load_snps_tbx().
#'
#' @param cnv see load_snps_tbx() documentation
#' @param samp see load_snps_tbx() documentation
#' @param snps see load_snps_tbx() documentation
#' @param in_out_ratio see load_snps_tbx() documentation
#' @param adjusted_lrr see load_snps_tbx() documentation
#' @param min_lrr see load_snps_tbx() documentation
#' @param max_lrr see load_snps_tbx() documentation
#' @param shrink_lrr see load_snps_tbx() documentation
#' @param w size of square side in pixel
#' @param z_ratio approximate ratio of z (the blank section between
#'        the two halves) compared to w
#' @param tmp_plot for developing, if set to 1 plot the "normal"
#'        LRR/BAF plot, if set to 2 plot the pixelated image in R
#'
#' @export
#'
#' @import data.table
#' @import ggplot2


plot_cnv <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                     w = 64, z_ratio = 0.1, tmp_plot = 0, min_lrr = -1.2, max_lrr = 1,
                     shrink_lrr = NULL) {
  # initial checks
  if ((w %% 2) != 0) stop('w must be even')

  # compute w, z and k values
  z <- round(w * z_ratio)
  if ((z %% 2) != 0) z <- z - 1
  k <- (w - z)/2
  if (((k*2) + z) != w) stop('something wrong')

  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # load snps data
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)
  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(data.table())
  }

  # move position to the x coordinates in the new system
  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)
  dt[, x := round(((position-ss)/(ee-ss)) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  dt_lrr <- copy(dt)
  dt_baf <- copy(dt)

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k)]
  dt_baf[, y := round(((baf-0)/(1-0)) * k) + k + z]

  # pixel coordinates must be > 0
  dt_baf[, ':=' (x = x+1, y = y+1)]
  dt_lrr[, ':=' (x = x+1, y = y+1)]
  w <- w+1

  # temporary check
  if (tmp_plot == 1) {
    a <- ggplot(dt_lrr, aes(x, y)) + geom_point() + xlim(0, w) + ylim(0, k) + theme_bw()
    b <- ggplot(dt_baf, aes(x, y)) + geom_point() + xlim(0, w) + ylim(k+z, w) + theme_bw()
    print(cowplot::plot_grid(b, a, ncol = 1))
  }

  # create the pixel map
  dt_lrr <- dt_lrr[, .N, by = c('x', 'y')]
  dt_baf <- dt_baf[, .N, by = c('x', 'y')]
  dt <- as.data.table(rbind(t(combn(1:w, 2)), t((combn(w:1, 2))))) # almost all combinations
  colnames(dt) <- c('x', 'y')
  dt <- rbind(dt, data.table(x = 1:w, y = 1:w)) # the diagonal was missing
  dt <- merge(dt, dt_lrr, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'a')
  dt <- merge(dt, dt_baf, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'b')
  dt[!is.na(a), value := a]
  dt[!is.na(b), value := b]
  dt[!is.na(a) & !is.na(b), value := a+b]
  dt[is.na(a) & is.na(b), value := 0]
  dt <- dt[, .(x, y, value)]

  if (tmp_plot == 2)
    print(ggplot(dt, aes(x, y, fill = value)) + geom_tile() + theme_bw() +
            scale_fill_gradient(low="white", high="black"))

  return(dt)

}
