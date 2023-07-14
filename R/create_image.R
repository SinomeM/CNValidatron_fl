#' Create a n*n machine-freindly image of the classic LRR BAF plots
#'
#' @param dt the output of `QCtreeCNV:::get_region_tabix()`
#' @param region one line data.table with CNV position information
#' @param n image dimension in pixels
#' @param n number of pixel per side of the final matrix
#' @param eps epsilon, number of empty pixels "between" the two plots, must be even
#' @param sides proportion of the CNV length in each side region (simmetric)
#' @param adj use `LRRadj` or not
#' @param path_img where to save the PNG
#'
#' @import data.table
#'
#' @export

cnv_image_matrix <- function(dt, region, n, eps = 4, sides = 0.5, adj = T, path_img) {

  # check if eps is even, n must be even as well
  if ((n %% 2) == 0 & (eps %% 2) != 0) stop("eps must be even if n is even")
  if ((n %% 2) != 0 & (eps %% 2) == 0) stop("eps must be odd if n is odd")

  if (eps > 0.1*n) stop("eps value should not exeed 10% of n")

  # extract the relevant markers
  ch <- region$chr; st <- region$start; en <- region$end; len <- en-st+1
  a <- st - (len*sides); b <- en + (len*sides)

  dt <- dt[chr == ch & between(position, a, b), ]

  if (adj) dt[, LRR := LRRadj]

  # LRR values higher than 2 are set to 2, lower than -2 to -2
  dt[LRR > 2, LRR := 2][LRR < -2, LRR := -2]

  lrr <- dt[ , .(position, LRR)]
  baf <- dt[ , .(position, BAF)]

  # two n*np tables descirbing the actual matrix, x is [0,n-1], y [0,np-1]
  np <- round((n-eps) / 2)
  lrr[, x := ((position-a) / (b-a)) * (n-1) ][, y := ((LRR+2) / 4) * (np-1)]
  baf[, x := ((position-a) / (b-a)) * (n-1) ][, y := BAF * (np-1)]

  npp <- n/2 + eps
  lrr[, y := y + npp]
  dtm <- round(rbind(lrr[, .(x, y)], baf[, .(x, y)]))

  # in the n*n table x and y are [1,n]
  dtm[, x := x+1][, y := y+1]

  # count the value of each "pixel" and scale it to [0,1]
  dtm <- dtm[, .N, by= list(x,y)]
  dtm[, N := ((N-min(dtm$N, na.rm = T)) / max(dtm$N, na.rm = T))]

  # create all possible x y positions, [1,n]
  tmp <- as.data.table(gtools::permutations(n, 2))
  colnames(tmp) <- c("x", "y")
  # the diagonal will be missing
  tmp <- rbind(tmp, data.table(x = 1:n, y = 1:n))

  # create the final table of pixel values and coordinates
  dtm <- merge(tmp, dtm, all.x=T)
  # change NAs to 0s
  dtm[is.na(N), N := 0]

  # from xy values to the actual matrix
  setorder(dtm, -y, x)
  dt <- matrix(dtm$N, ncol = n, byrow = T)

  png::writePNG(dt, path_img)

  return(dt)
}
