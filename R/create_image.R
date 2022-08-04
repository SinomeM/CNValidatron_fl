#' Create a n*n machine-freindly image of the classic LRR BAF plots
#'
#' @param n image dimension in pixels
#' @dt the output of `QCtreeCNV:::get_region_tabix()`
#' @region one line data.table with CNV position information
#' @n number of pixel per side of the final matrix
#' @eps epsilon, number of empty pixels "between" the two plots
#' @sides proportion of the CNV length in each side region (simmetric)
#'
#' @import data.table
#'
#' @export

cnv_image_matrix <- function(dt, region, n, eps = 4, sides = 0.5) {

  # extract the relevant markers
  ch <- region$chr; st <- region$start; en <- region$end; len <- en-st
  a <- st - (len*sides); b <- en + (len*sides)

  dt <- dt[chr == ch & between(position, a, b), ]

  lrr <- dt[ , .(position, LRR)]
  baf <- dt[ , .(position, BAF)]

  # two n*np matrix
  np <- (n/2) - (eps/2)

  lrr[, x := ((position-a) / (b-a)) * n ][, y := ((LRR+2) / 4) * np]
  baf[, x := ((position-a) / (b-a)) * n ][, y := BAF * np]

  # compose the final matrix
  npp <- n/2 + eps
  lrr[, y := y + npp]

  dtm <- round(rbind(lrr[, .(x, y)], baf[, .(x, y)]))

  dtm <- dtm[, .N, by= list(x,y)]

  # scale it to 0:128
  dtm[, N := round(((N-min(dtm$N)) / max(dtm$N)) * 128)]

  ggplot(dtm, aes(x,y, fill = -N)) +
    geom_tile(show.legend=F) +
    #geom_tile() +
    theme(axis.title = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), panel.background = element_blank()) +
    ylim(0, n) + xlim(0,n) +
    scale_fill_distiller(type = "seq", palette = "Greys")

}
