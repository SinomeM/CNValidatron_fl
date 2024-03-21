#' Load intensity data in R from tabix file
#'
#' A function to lead the snps data from the tabix indexed intensity file.
#' It will also perform some preprocessing, e.g. reduce the LRR interval to [-1.4 , 1.2],
#'
#' @param cnv one line data.table in the usual cnv format
#' @param samp one line samples file in the usual format
#' @param snps the snps file for PennCNV, as data.table
#' @param in_out_ratio ratio of bp outside the cnv vs inside, per side. A value of one means
#'        the CNV length will be added on each side
#' @param adjusted_lrr load GC-adjusted LRR on normal LRR (T/F)
#' @param min_lrr minimum LRR value (lower values are set to the minimum)
#' @param max_lrr maximum LRR value
#' @param shrink_lrr shrink LRR values toword the mean (done separately for
#'        SNPs before / in / after the CNV
#'
#' @export
#'
#' @import data.table


load_snps_tbx <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                          min_lrr = -1.2, max_lrr = 1, shrink_lrr = NULL) {
  chr <- cnv$chr
  start <- cnv$start
  end <- cnv$end
  len <- end - start + 1
  tbx_path <- samp$file_path_tabix

  st <- start - (in_out_ratio*len)
  st <- ifelse(st < 0, 0, st)

  dt <- fread(cmd = paste0("tabix ", tbx_path, " ", chr, ":", st,
                          "-", end + (in_out_ratio*len)), header = F)

  if (nrow(dt) == 0) {
    warning('File: ', tbx_path, ' seems empty or broken\n')
    return(data.table())
  }

  if (ncol(dt) == 7) colnames(dt) <- c("chr", "position", "end", "LRR", "LRRadj", "BAF", "snp")
  else {
    if (adjusted_lrr) colnames(dt) <- c("chr", "position", "end", "LRR", "BAF", "LRRadj")
    else colnames(dt) <- c("chr", "position", "end", "LRR", "BAF")
  }

  if (!is.null(snps)) dt <- dt[paste0(chr, position) %in% snps[, paste0(Chr, Position)], ]

  setorder(dt, position)

  if (adjusted_lrr) setnames(dt, c('BAF', 'LRRadj'), c('baf', 'lrr'))
  else setnames(dt, c('BAF', 'LRR'), c('baf', 'lrr'))

  ## some preprocessing ##
  # restrict the lrr space
  dt[lrr > max_lrr, lrr := max_lrr][lrr < min_lrr, lrr := min_lrr]
  # id lrr is missing outside of the cnv 'impute' it
  dt[is.na(lrr) & !between(position, start, end), lrr := dt[!between(position, start, end),
                                                              mean(lrr, na.rm = T)]]
  # if baf is missing outside the cnv set it to either 0 or 1
  dt[is.na(baf) & !between(position, start, end), baf := sample(rep(0:1, length.out = .N))]
  # if lrr or baf is missing inside the cnv exclude the point
  dt <- dt[!((is.na(lrr) | is.na(baf)) & between(position, start, end)),]

  ## some more processing ##
  dt[position < start, group := 1][between(position, start, end), group := 2][
       position > end, group := 3]
  ms1 <- dt[group == 1, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms2 <- dt[group == 2, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms3 <- dt[group == 3, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]

  if (!is.null(shrink_lrr)) {
    # snps in each group get pulled towards the group mean relative to their distance
    dt[group == 1 & lrr > ms1[1], lrr := lrr - (lrr-ms1[1])*shrink_lrr][
         group == 1 & lrr < ms1[1], lrr := lrr + (lrr-ms1[1])*shrink_lrr]
    dt[group == 2 & lrr > ms2[1], lrr := lrr - (lrr-ms2[1])*shrink_lrr][
         group == 2 & lrr < ms2[1], lrr := lrr + (lrr-ms2[1])*shrink_lrr]
    dt[group == 3 & lrr > ms3[1], lrr := lrr - (lrr-ms3[1])*shrink_lrr][
         group == 3 & lrr < ms3[1], lrr := lrr + (lrr-ms3[1])*shrink_lrr]
  }

  # ouliers removal, 3SDs
  dt[group == 1 & !between(lrr, ms1[1]-3*ms1[2], ms1[1]+3*ms1[2]), lrr := NA]
  dt[group == 1 & !between(lrr, ms2[1]-3*ms2[2], ms2[1]+3*ms2[2]), lrr := NA]
  dt[group == 1 & !between(lrr, ms3[1]-3*ms3[2], ms3[1]+3*ms3[2]), lrr := NA]
  dt <- dt[!is.na(lrr), ]



  return(dt[, .(chr, position, lrr, baf)])

}
