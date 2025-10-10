#' Create the matrix for the PNG image
#'
#' This function create the pixel matrix that can be saved
#' as a PNG for further use. It will also perform some preprocessing, e.g.
#' reduce the LRR interval to [-1.4 , 1.2],
#'
#' @param cnv see load_snps_tbx() documentation
#' @param samp see load_snps_tbx() documentation
#' @param shrink_lrr shrink LRR values toword the mean (done separately for
#'        SNPs before / in / after the CNV
#' @param min_lrr minimum LRR value (lower values are set to the min_lrr)
#' @param max_lrr maximum LRR value (higher values are set to max_lrr)
#' @param tmp_plot for developing, if set to 1 plot the "normal"
#'        LRR/BAF plot, if set to 2 plot the pixelated image in R
#'
#' @export
#'
#' @import data.table
#' @import ggplot2


plot_cnv <- function(cnv, samp, tmp_plot = 0, shrink_lrr = 0.2, tabix_data = NULL,
                     min_lrr = -1.4, max_lrr = 1.3, blank_small = F,
                     # the following parameters should not be changed by most users
                     w = 96, z = 4, k1 = 31, k2 = 26,
                     l_wind = 20000000, # top row Mbp
                     mx_lr = 2) {  # top row LRR range

  dt <- tabix_data
  # w k1, k2, z and in_out_ratio are fixed for the moment

  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # in_out_ratio is now a function of the legnth
  len <- cnv$end - cnv$start + 1
  if (len <= 100000) in_out_ratio <- 9
  if (between(len, 100001, 1000000)) in_out_ratio <- 7
  if (len > 1000000) in_out_ratio <- 5
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)

# from load_snps_tbx() .......
  start <- cnv$start
  end <- cnv$end
  len <- end - start + 1

  st <- start - (in_out_ratio*len)
  st <- ifelse(st < 0, 0, st)

  dt_whole <- copy(dt)
  dt <- dt[between(position, start - len*in_out_ratio, end + len*in_out_ratio), ]

  ## some preprocessing ##
  # restrict the lrr space
  dt[lrr > max_lrr, lrr := max_lrr][lrr < min_lrr, lrr := min_lrr]
  # if lrr or baf is missing exclude the point
  dt <- dt[!(is.na(lrr) | is.na(baf)),]

  # compute mean and SD in these three groups, could be simplified using dt[,,by]
  dt[position < start, group := 1][
     between(position, start, end), group := 2][position > end, group := 3]
  ms1 <- dt[group == 1, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms2 <- dt[group == 2, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms3 <- dt[group == 3, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]

  if (!is.null(shrink_lrr)) {
    # snps in each group get pulled towards the group mean proportionally
    # to their distance and shrink_lrr
    dt[group == 1 & lrr > ms1[1], lrr := lrr - abs(lrr-ms1[1])*shrink_lrr][
         group == 1 & lrr < ms1[1], lrr := lrr + abs(ms1[1]-lrr)*shrink_lrr]

    dt[group == 2 & lrr > ms2[1], lrr := lrr - abs(lrr-ms2[1])*shrink_lrr][
         group == 2 & lrr < ms2[1], lrr := lrr + abs(ms2[1]-lrr)*shrink_lrr]

    dt[group == 3 & lrr > ms3[1], lrr := lrr - abs(lrr-ms3[1])*shrink_lrr][
         group == 3 & lrr < ms3[1], lrr := lrr + abs(ms3[1]-lrr)*shrink_lrr]
  }

  # ouliers removal, 3SDs. Should not do anything after the rest
  dt[group == 1 & !between(lrr, ms1[1]-3*ms1[2], ms1[1]+3*ms1[2]), lrr := NA]
  dt[group == 2 & !between(lrr, ms2[1]-3*ms2[2], ms2[1]+3*ms2[2]), lrr := NA]
  dt[group == 3 & !between(lrr, ms3[1]-3*ms3[2], ms3[1]+3*ms3[2]), lrr := NA]
  dt <- dt[!is.na(lrr), ]
# ........


  if (nrow(dt) == 0) {
    warning('Empty tabix, no image generated for sample', samp$sample_ID)
    return(list(data.table(), 0))
  }

  n_real_snps <- dt[between(position, cnv$start, cnv$end), .N]

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

  dt_whole <- dt_whole[between(position, ss2, ee2),]
  dt_whole[, x := round(((position-ss2)/(ee2-ss2)) * w)]
  dt_whole[, y := round(((lrr-(-mx_lr))/(mx_lr-(-mx_lr))) * k2) + (k1*2 + z*2)]
  dt_whole[, ':=' (x = x+1, y = y+1)]

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
    c <- ggplot(dt_whole, aes(position, lrr)) + geom_point(alpha = 0.1, colour = 'purple') +
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
    c <- ggplot(dt_whole, aes(x, y)) + geom_point(alpha = 0.05, colour = 'purple') +
           xlim(0, w) + ylim((k1*2)+(z*2) + 1, w + 1) + theme_bw() +
           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(), axis.title.y = element_blank())
    return(cowplot::plot_grid(c, b, a, ncol = 1))
  }


  # create the pixel map
    dt_lrr <- get_normalised_pixel_values(dt_lrr, w, in_out_ratio)
    dt_baf <- get_normalised_pixel_values(dt_baf, w, in_out_ratio)
    dt_whole <- get_normalised_pixel_values(dt_whole, w, in_out_ratio)

  if (blank_small) {
    dt_bf <- dt_lrr[x %in% 47:50, ]
    # here we are already in normalised pixel intensities
    if (dt_lrr[N >= 0, .N] < some_number)
      message('there are very few point available for the CNVs')
    # consider exuding the call, e.g. make the PNG blank and train the model
    # to call blank images as "too small"
  }


  # create the full picture
  dt <- as.data.table(rbind(t(combn(1:w, 2)), t((combn(w:1, 2))))) # almost all combinations
  colnames(dt) <- c('x', 'y')
  dt <- rbind(dt, data.table(x = 1:w, y = 1:w)) # the diagonal was missing
  # fill in the values
  dt <- merge(dt, dt_lrr, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'a')
  dt <- merge(dt, dt_baf, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'b')
  dt <- merge(dt, dt_whole, by = c('x', 'y'), all.x = T)
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


  return(list(dt, n_real_snps))

}


get_normalised_pixel_values <- function(dt, w, in_out_ratio) {

  if (in_out_ratio == 9) steps <- 2
  if (in_out_ratio == 7) steps <- 3
  if (in_out_ratio == 5) steps <- 5

  wind_size <- w %/% 2^steps

  # get N for each pixel
  dt <- dt[, .N, by = c('x', 'y')]
  dt[N >= 1, N := N + 1]
  dt[, N := log(N + 0.0000000001)]

  for (i in 1:steps-1) {
    # create windows
    wind_size <- wind_size * i
    dt[, wind := x %/% wind_size]
    # z score normalisation per window
    dt[, N := (N - mean(N)) / sd(N), by = wind]
  }
  # min max normalisation, move to [0,1]. The minimum is moved to 0.05 to increase the contrast at low N
  dt[, N := (((N-min(N)) / (max(N)-min(N))) / 1.081) + 0.075, by = wind]

  return(dt)
}

get_normalised_pixel_values_simple <- function(dt) {

  # get N for each pixel
  dt <- dt[, .N, by = c('x', 'y')]
  # increase the contrast between "small N" and 0
  dt[N >= 1, N := N + 1]

  # min max normalization on the log scale
  dt[, N := log(N + 0.0000000001)]
  dt[, N := (((N-min(N)) / (max(N)-min(N))) / 1.081) + 0.075, ]

  return(dt)
}
