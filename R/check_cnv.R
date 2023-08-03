
check_cnv <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                      min_lrr = -1.2, max_lrr = 1, shrink_lrr = NULL) {

  # load snps data
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)

  dt[between(position, cnv$start, cnv$end), inside := T][is.na(inside), inside := F]

  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len)
  ee <- cnv$end + (in_out_ratio*len)

  a <- ggplot(dt, aes(x = position, y = baf, colour = inside)) + geom_point() +
    xlim(ss, ee) + theme_bw()
  b <- ggplot(dt, aes(x = position, y = lrr, colour = inside)) + geom_point() +
    xlim(ss, ee) + ylim(min_lrr, max_lrr) + theme_bw()

  cowplot::plot_grid(a,b, ncol = 1)

}

