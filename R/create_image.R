#' Create a n*n machine-freindly image of the classic LRR BAF plots
#'
#' @param n image dimension in pixels
#' @dt the output of `QCtreeCNV:::get_region_tabix()`
#' @region one line data.table with CNV position information
#' @n number of pixel per side of the final matrix
#' @eps epsilon, number of empty pixels "between" the two plots, must be even
#' @sides proportion of the CNV length in each side region (simmetric)
#'
#' @import data.table
#'
#' @export

cnv_image_matrix <- function(dt, region, n, eps = 4, sides = 0.5, adj = T, path_img) {

  # check if eps is even, n must be even as well
  if ((n %% 2) == 0) & (eps %% 2) != 0) stop("eps must be even if n is even")
  if ((n %% 2) != 0) & (eps %% 2) == 0) stop("eps must be odd if n is odd")

  # extract the relevant markers
  ch <- region$chr; st <- region$start; en <- region$end; len <- en-st+1
  a <- st - (len*sides); b <- en + (len*sides)

  dt <- dt[chr == ch & between(position, a, b), ]
  
  if (adj) dt[, LRR := LRRadj]

  # LRR values higher than 2 are set to 2, lower than -2 to -2
  dt[LRR > 2, LRR := 2][LRR < -2, LRR := -2]

  lrr <- dt[ , .(position, LRR)]
  baf <- dt[ , .(position, BAF)]

  # two n*np matrix
  np <- round((n-eps) / 2)

  lrr[, x := ((position-a) / (b-a)) * n ][, y := ((LRR+2) / 4) * np]
  baf[, x := ((position-a) / (b-a)) * n ][, y := BAF * np]

  # compose the final matrix
  npp <- n/2 + eps
  lrr[, y := y + npp]

  dtm <- round(rbind(lrr[, .(x, y)], baf[, .(x, y)]))

  dtm <- dtm[, .N, by= list(x,y)]

  # scale it to [0,1]
  dtm[, N := ((N-min(dtm$N, na.rm = T)) / max(dtm$N, na.rm = T))]

  # create all possible x y positions
  tmp <- as.data.table(gtools::permutations(n, 2))
  colnames(tmp) <- c("x", "y")
  # the diagonal will be missing
  tmp <- rbind(tmp, data.table(x = 1:n, y = 1:n))

  dtm <- merge(tmp, dtm, all.x=T)
  # change NAs to 0s
  dtm[is.na(N), N := 0]

  # from xy values to the actual matrix
  setorder(dtm, -y, x)
  dt <- matrix(dtm$N, ncol = n, byrow = T)

  png::writePNG(dt, path_img)

  return(dt)
}
