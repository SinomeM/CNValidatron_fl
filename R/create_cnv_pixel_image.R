#' Create the matrix for the PNG image
#'
#' This function create the pixel matrix that can be saved
#' as a PNG for further use. All processing is done in the
#' function load_snps_tbx().
#'
#' @param cnv see load_snps_tbx() documentation
#' @param samp see load_snps_tbx() documentation
#' @param snps see load_snps_tbx() documentation
#' @param adjusted_lrr see load_snps_tbx() documentation
#' @param min_lrr see load_snps_tbx() documentation
#' @param max_lrr see load_snps_tbx() documentation
#' @param shrink_lrr see load_snps_tbx() documentation
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
  w <- 96
  z <- 4
  k1 <- 31
  k2 <- 26
  in_out_ratio <- 5
  # top row Mbp
  l_wind <- 20000000
  # top row LRR range
  mx_lr = 2.5
  # top row normalisation windows size
  norm_wind_size = 6

  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # load snps data, ALL chromosome is loaded now!
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)
  # keep the full chromsome for the third row of the png
  dt_big <- dt[[2]]
  dt <- dt[[1]]
  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(data.table())
  }

  # bottom and middle row, select the relevant points and
  # move position to the x coordinates in the new system
  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)
  dt <- dt[between(position, ss, ee), ]
  dt[, x := round(((position-ss)/(ee-ss)) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  # formally the second copy is not needed
  dt_lrr <- copy(dt)
  dt_baf <- copy(dt)

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k1)]
  dt_baf[, y := round(((baf-0)/(1-0)) * k1) + k1 + z]

  # pixel coordinates must be > 0
  dt_baf[, ':=' (x = x+1, y = y+1)]
  dt_lrr[, ':=' (x = x+1, y = y+1)]

  # top row
  center <- cnv$start + len/2
  a <- center - l_wind/2
  b <- center + l_wind/2
  # if a is out of bound then the 30Mbp region must start from the first SNP
  if (a < dt_big[, min(position)]) {
    a <- dt_big[, min(position)]
    b <- a+l_wind
  }
  else
  # same thing on the other side, both cannot be true at the same time in a human genome
    if (b > dt_big[, max(position)]) {
      b <- dt_big[, max(position)]
      a <- b-l_wind
    }

  dt_big <- dt_big[between(position, a, b),]
  dt_big[, x := round(((position-a)/(b-a)) * w)]
  dt_big[, y := round(((lrr-(-mx_lr))/(mx_lr-(-mx_lr))) * k2) + (k1*2 + z*2)]
  dt_big[, ':=' (x = x+1, y = y+1)]

  w <- w+1

  # plot in the original space
  if (tmp_plot == 1) {
    a <- ggplot(dt_lrr, aes(position, lrr)) + geom_point(alpha = 0.3, colour = 'red') +
           ylim(min_lrr, max_lrr) + theme_bw() +
           geom_segment(x = cnv$start, xend = cnv$end, y = 0, yend = 0, linetype = 3) +
           theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    b <- ggplot(dt_baf, aes(position, baf)) + geom_point(alpha = 0.3, colour = 'blue') +
           theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    c <- ggplot(dt_big, aes(position, lrr)) + geom_point(alpha = 0.1, colour = 'purple') +
           geom_segment(x = cnv$start, xend = cnv$end, y = 0, yend = 0) +
           ylim(min_lrr, max_lrr) + theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    return(cowplot::plot_grid(c, b, a, ncol = 1))
  }

  # plot in the new coordinates but suitable for humans
  if (tmp_plot == 2) {
    a <- ggplot(dt_lrr, aes(x, y)) + geom_point(alpha = 0.1, colour = 'red') +
           xlim(0, w) + ylim(0+1, k1+1) + theme_bw() +
           theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    b <- ggplot(dt_baf, aes(x, y)) + geom_point(alpha = 0.1, colour = 'blue') +
           xlim(0, w) + ylim(k1+z + 1, k1*2 + z + 1) + theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    c <- ggplot(dt_big, aes(x, y)) + geom_point(alpha = 0.05, colour = 'purple') +
           xlim(0, w) + ylim((k1*2)+(z*2) + 1, w + 1) + theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    return(cowplot::plot_grid(c, b, a, ncol = 1))
  }

  # create the pixel map, and normalise values between 0 and 1
  dt_lrr <- dt_lrr[, .N, by = c('x', 'y')][, N := N/max(N)]
  dt_baf <- dt_baf[, .N, by = c('x', 'y')][, N := N/max(N)]
  # the top row normalisation is performed in windows first, then on all
  dt_big <- dt_big[, .N, by = c('x', 'y')]
  dt_big[, wind1 := x %% norm_wind_size]
  dt_big[, wind2 := x %% (norm_wind_size+2)]
  dt_big[, n := N/max(N), by = wind1]
  dt_big[, n := N/max(N), by = wind2]
  dt_big[, N := n][, n := NULL][, N := N/max(N)]

  # create the full picture
  dt <- as.data.table(rbind(t(combn(1:w, 2)), t((combn(w:1, 2))))) # almost all combinations
  colnames(dt) <- c('x', 'y')
  dt <- rbind(dt, data.table(x = 1:w, y = 1:w)) # the diagonal was missing
  # fill in the values
  dt <- merge(dt, dt_lrr, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'a')
  dt <- merge(dt, dt_baf, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'b')
  dt <- merge(dt, dt_big, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'c')
  dt[!is.na(a), value := a]
  dt[!is.na(b), value := b]
  dt[!is.na(c), value := c]
  # dt[!is.na(a) & !is.na(b), value := a+b] # this should not be a case
  dt[is.na(a) & is.na(b) & is.na(c), value := 0]
  dt <- dt[, .(x, y, value)]

  # plot almost identical to the final PNG
  if (tmp_plot == 3)
    return(ggplot(dt, aes(x, y, fill = value)) + geom_tile() + theme_minimal() +
             scale_fill_gradient(low="white", high="black"))

  dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis

  return(dt)

}
