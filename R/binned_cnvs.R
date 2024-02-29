#' Compute the overlap between a CNV table and a list of genomic bins
#'
#' This can be regular windows throughout the genome or exons / genes for example.
#'
#' @param cnvs usual CNV `data.table`
#' @param format format for the output table, "wide" and "biot" return a list
#' @param bins for computing the overlap, by default are regular bins created by the helper function `binned_genome`
#'
#' @import data.table
#' @export

binned_cnvs <- function(cnvs, format = c('long', 'count', 'wide', 'both'),
                        bins = binned_genome(chrs_sten = QCtreeCNV::hg19_start_end_centromeres, bin_size = 100000, chrs = 1:22)) {

  dt <- bins

  # compute the table in "long format", one line per CNV-bin match
  dto <- data.table()
  for (i in dt[, unique(chr)]) {
    dtx <- cnvs[chr == i, .(chr, start, end, sample_ID, GT)]
    dty <- dt[chr == i, ]
    setkey(dtx, start, end)
    setkey(dty, start, end)

    dto <- rbind(dto, unique(foverlaps(dtx, dty)))
  }

  if (format == 'long') return(dto)
  # add the count table, one line per "variant"
  if (format == 'count') return(list(dto, dto[, .N, by = c('chr', 'start', 'end', 'GT')]))

  # create the table in the wide format, one column per sample
  dtow <- dcast(dto, ix + chr + start + end ~ sample_ID, value.var = 'GT')
  gc()
  if (format == 'wide') return(dtow)
  if (format == 'both') return(list(dto, dtow))
}

binned_genome <- function(chrs_sten = QCtreeCNV::hg19_start_end_centromeres,
                          bin_size = 500000, chrs = 1:22) {
  # create the binned genome
  dt <- data.table()
  for (i in chrs_sten[chr %in% chrs, unique(chr)]) {
    lims <- chrs_sten[chr == i, c(start, end)]
    n <- round((lims[2] - lims[1] + 1) / bin_size) + 1
    tmp <- data.table(ix = 1:n, chr = as.integer(i))
    tmp[, ':=' (start = 0 + (bin_size * (ix-1)), end = bin_size - 1 + (bin_size * (ix-1)))]
    tmp[end > lims[2], end := lims[2]]
    dt <- rbind(dt, tmp)
  }
  return(dt)
}
