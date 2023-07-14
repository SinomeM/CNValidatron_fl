
# A function to lead the snps data from the tabix indexed intensity file.
# It will also perform some preporcessing, e.g. reduce the LRR interval to [-1.4 , 1.2],
# in hte future some denoising may be added

# @param cnv one line data.table in the usual cnv format
# @param samp one line samples file in the usual format
# @param snps the snps file for PennCNV, as data.table
# @param in_out_ratio ratio of bp otside the cnv vs inside, per side. A value of one means
#        the cnv length will be addded on each side
load_snps_tbx <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T) {
  chr <- cnv$chr
  start <- cnv$start
  end <- cnv$end
  len <- end - start + 1
  tbx_path <- samp$file_path_tabix

  dt <- fread(cmd = paste0("tabix ", tbx_path, " ", chr, ":", start - (in_out_ratio*len),
                          "-", end + (in_out_ratio*len)), header = F)

  if (adjusted_lrr) {
    colnames(dt) <- c("chr", "position", "end", "LRR", "BAF", "LRRadj")
    dt[, end := NULL]
  }
  else colnames(dt) <- c("chr", "position", "LRR", "BAF")

  if (!is.null(snps)) dt <- dt[paste0(chr, position) %in% snps[, paste0(Chr, Position)], ]

  setorder(dt, position)

  if (adjusted_lrr) setnames(dt, c('BAF', 'LRRadj'), c('baf', 'lrr'))
  else setnames(dt, c('BAF', 'LRR'), c('baf', 'lrr'))

  # some preprocessing
  # restrict the lrr space
  dt[lrr > 1.2, lrr := 1.2][lrr < -1.4, lrr := -1.4]
  # id lrr is missing outside of the cnv 'impute' it
  dt[is.na(lrr) & !between(position, start, end), lrr := dt[!between(position, start, end), mean(lrr)]]
  # if baf is missing outised the cnv set it to either 0 or 1
  dt[is.na(baf) & !between(position, start, end), sample(rep(0:1, length.out = .N))]
  # if lrr or baf is missing inside the cnv exclude the point
  dt <- dt[!((is.na(lrr) | is.na(baf)) & between(position, start, end)),]

  return(dt)

}

plot_cnv <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                     w = 128, z_ratio = 0.1) {
  # initial checks
  if ((w %% 2) != 0) stop('w must be even')

  # compute w, z and k values
  z <- round(w * z_ratio)
  if ((z %% 2) != 0) z <- z - 1
  k <- (w - z)/2
  if ((k %% 2) != 0) stop('something wrong')
  if (((k*2) + z) != w) stop('something else wrong')

  # load snps data
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr)

  # move position to the x coordinates in the new system
  dt[, x := round(((position-min(position))/(max(position)-min(position))) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  dt_lrr <- copy(dt)
  dt_baf <- copy(dt)

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-min(lrr))/(max(lrr)-min(lrr))) * k)]
  dt_baf[, y := round(((baf-min(baf))/(max(baf)-min(baf))) * k) + k + z]

  return(list(dt_lrr, dt_baf))

}
