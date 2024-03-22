#' Plot LRR/BAF for a given CNV
#'
#' This function is mostly needed for developing to check how CNVs looks like
#' "normally" and compare with the PNG form.
#'
#' @param cnv see load_snps_tbx() documentation
#' @param samp see load_snps_tbx() documentation
#' @param snps see load_snps_tbx() documentation
#' @param in_out_ratio see load_snps_tbx() documentation
#' @param adjusted_lrr see load_snps_tbx() documentation
#' @param min_lrr see load_snps_tbx() documentation
#' @param max_lrr see load_snps_tbx() documentation
#' @param shrink_lrr see load_snps_tbx() documentation
#'
#' @import data.table
#' @import ggplot2

check_cnv <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                      min_lrr = -1.2, max_lrr = 1, shrink_lrr = NULL) {

  # load snps data, check behaviour if needed. load_snps_tbx() has changed
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)[[1]]

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

