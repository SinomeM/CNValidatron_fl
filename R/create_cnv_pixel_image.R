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
                     # the following parameters should not be changed by most users
                     shrink_lrr = 0.1, w = 96, z = 4, k1 = 31, k2 = 26,
                     l_wind = 20000000, # top row Mbp
                     mx_lr = 2) {  # top row LRR range

  # w k1, k2, z and in_out_ratio are fixed for the moment


  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # in_out_ratio is now a function of the legnth
  len <- cnv$end - cnv$start + 1
  if (len <= 100000) in_out_ratio <- 9
  if (between(len, 100001, 1000000)) in_out_ratio <- 7
  if (len > 1000000) in_out_ratio <- 5
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)


  # load snps data, ALL chromosome is loaded now!
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)
  # keep the full chromsome for the third row of the png
  dt_big <- dt[[2]]
  dt <- dt[[1]]
  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(data.table())
  }


  # bottom and middle row, move position to the x coordinates in the new system
  # dt <- dt[between(position, ss, ee), ] # already done in load_snps_tbx()
  dt[, x := round(((position-ss)/(ee-ss)) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  dt_lrr <- copy(dt)
  dt_baf <- dt

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k1)]
  dt_baf[, y := round(((baf-0)/(1-0)) * k1) + k1 + z]

  # pixel coordinates must be > 0
  dt_baf[, ':=' (x = x+1, y = y+1)]
  dt_lrr[, ':=' (x = x+1, y = y+1)]


  # top row
  #  the top row is zoomed out by a factor of 3
  ss2 <- cnv$start - (in_out_ratio*len*3);  ee2 <- cnv$end + (in_out_ratio*len*3)
  # if it is not at least 12.5Mbp then force it to 12.5Mbp
  len_diff <- (12500000 - (ee2 - ss2 + 1)) / 2
  if (len_diff > 0) {
    ss2 <- ss2 - len_diff
    ee2 <- ee2 + len_diff
  }

  dt_big <- dt_big[between(position, ss2, ee2),]
  dt_big[, x := round(((position-ss2)/(ee2-ss2)) * w)]
  dt_big[, y := round(((lrr-(-mx_lr))/(mx_lr-(-mx_lr))) * k2) + (k1*2 + z*2)]
  dt_big[, ':=' (x = x+1, y = y+1)]

  w <- w+1

  # plot in the original space
  if (tmp_plot %in% c(1,3)) {
    a <- ggplot(dt_lrr, aes(position, lrr)) + geom_point(alpha = 0.3, colour = 'red') +
           ylim(min_lrr, max_lrr) + theme_bw() + xlim(ss, ee) +
           geom_segment(x = cnv$start, xend = cnv$end, y = 0, yend = 0, linetype = 3) +
           theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    b <- ggplot(dt_baf, aes(position, baf)) + geom_point(alpha = 0.3, colour = 'blue') +
           theme_bw() + xlim(ss, ee) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    c <- ggplot(dt_big, aes(position, lrr)) + geom_point(alpha = 0.1, colour = 'purple') +
           geom_segment(x = cnv$start, xend = cnv$end, y = 0, yend = 0) +
           ylim(min_lrr, max_lrr) + theme_bw() + xlim(ss2, ee2) +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())

    pl <- cowplot::plot_grid(c, b, a, ncol = 1)
    if (tmp_plot == 1) return(pl)
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


  # create the pixel map
  dt_lrr <- get_normalised_pixel_values(dt_lrr, w, in_out_ratio)
  dt_baf <- get_normalised_pixel_values(dt_baf, w, in_out_ratio)
  dt_big <- get_normalised_pixel_values(dt_big, w, in_out_ratio)


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
  if (tmp_plot == 3) {
    pl2 <- ggplot(dt, aes(x, y, fill = value)) + geom_tile() + theme_minimal() +
             scale_fill_gradient(low="white", high="black")
    return(cowplot::plot_grid(pl, pl2, ncol = 2))
  }

  dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis

  return(dt)

}


get_normalised_pixel_values <- function(dt, w, in_out_ratio) {

  if (in_out_ratio == 9) steps <- 2
  if (in_out_ratio == 7) steps <- 3
  if (in_out_ratio == 5) steps <- 5

  wind_size <- w %/% 2^steps

  # get N for each pixel
  dt <- dt[, .N, by = c('x', 'y')]
  dt[, N := as.numeric(N)]

  for (i in 1:steps-1) {
    # create windows
    wind_size <- wind_size * i
    dt[, wind := x %/% wind_size]
    # z score normalisation per window
    dt[, N := (N - mean(N)) / sd(N), by = wind]
  }
  # min max normalisation, move to [0,1]. The minimum is moved to 0.1 to increase the contrast at low N
  dt[, N := (((N-min(N)) / (max(N)-min(N))) / 1.11) + 0.1, by = wind]

  return(dt)
}
