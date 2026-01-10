#' Load intensity data in R from tabix file
#'
#' A function to lead the snps data from the tabix indexed intensity file.
#'
#' @param cnv one line data.table in the usual cnv format
#' @param samp one line samples file in the usual format
#' @param snps the snps file for PennCNV, as data.table
#' @param adjusted_lrr load GC-adjusted LRR on normal LRR (T/F)
#' @param shrink_lrr 
#'
#' @export
#'
#' @import data.table


load_snps_tbx <- function(cnv, samp, snps = NULL, adjusted_lrr = T) {
  chr <- cnv$chr
  tbx_path <- samp$file_path_tabix

  # load the whole chromosome now
  dt <- fread(cmd = paste0("tabix ", tbx_path, " ", chr), header = F)

  if (nrow(dt) == 0) {
    warning('File: ', tbx_path, ' seems empty or broken\n')
    return(data.table())
  }

  if (adjusted_lrr) colnames(dt) <- c("chr", "position", "end", "LRR", "BAF", "LRRadj")
  else colnames(dt) <- c("chr", "position", "end", "LRR", "BAF")

    # filter SNPs if snp object is provided
  if (!is.null(snps))
    dt <- dt[paste0(chr, position) %in% snps[, paste0(Chr, Position)], ]

  setorder(dt, position)

  if (adjusted_lrr) setnames(dt, c('BAF', 'LRRadj'), c('baf', 'lrr'))
  else setnames(dt, c('BAF', 'LRR'), c('baf', 'lrr'))

  return(dt[, .(chr, position, lrr, baf)])
}
