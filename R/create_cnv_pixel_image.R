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


plot_cnv <- function(cnv, samp, snps = NULL, adjusted_lrr = T,
                     tmp_plot = 0, min_lrr = -1.4, max_lrr = 1.3,
                     shrink_lrr = NULL) {

  # w k1, k2, z and in_out_ratio are fixed for the moment
  w <- 100
  z <- 4
  k1 <- 32
  k2 <- 28
  in_out_ratio <- 4

  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # load snps data, ALL chromosome is loaded now!
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)
  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(data.table())
  }

  # keep the full chromsome for the third row of the png
  dt_big <- copy(dt)

  # bottom and middle row, select the relevant points and
  # move position to the x coordinates in the new system
  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)
  dt <- dt[between(position, ss, en), ]
  dt[, x := round(((position-ss)/(ee-ss)) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  dt_lrr <- copy(dt)
  dt_baf <- copy(dt)

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k1)]
  dt_baf[, y := round(((baf-0)/(1-0)) * k1) + k1 + z]

  # pixel coordinates must be > 0
  dt_baf[, ':=' (x = x+1, y = y+1)]
  dt_lrr[, ':=' (x = x+1, y = y+1)]

  # top row, ~30 Mbp LRR
  center <- cnv$start + len/2
  a <- center - 15000000
  b <- center + 15000000
  # if a is out of bound then the 30Mbp region must start from the first SNP
  if (a < dt_big[, min(position)]) {
    a <- dt_big[, min(position)]
    b <- a+30000000
  }
  else
  # same thing on the other side, both cannot be true at the same time in a human genome
    if (b > dt_big[, max(position)]) {
      b <- dt_big[, max(position)]
      a <- b-30000000
    }

  dt_big <- dt_big[between(position, a, b),]
  dt_big[, x := round(((position-a)/(b-a)) * w)]
  dt_big[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k2) + (k1*2 + z*2)]
  dt_big[, ':=' (x = x+1, y = y+1)]
  ## HERE ##


  w <- w+1

  # temporary check
  if (tmp_plot == 1) {
    a <- ggplot(dt_lrr, aes(x, y)) + geom_point() + xlim(0, w) + ylim(0, k) + theme_bw()
    b <- ggplot(dt_baf, aes(x, y)) + geom_point() + xlim(0, w) + ylim(k1+z, k1*2 + z) + theme_bw()
    c <- ggplot(dt_big, aes(x, y)) + geom_point() + xlim(0, w) + ylim((k1*2)+(z*2), w) + theme_bw()
    return(cowplot::plot_grid(b, a, c, ncol = 1))
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

  dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis

  return(dt)

}
