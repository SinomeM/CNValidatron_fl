#' Quick fisher test association function
#'
#' Simple function to perform a fisher test for each
#' marker against a phenotype.
#'
#'
#'
#'
#' @import data.table


quick_assoc <- function(markers, pheno, test = 'fisher') {
  # check minimal columns
  if (!all(c('sample_ID', 'marker_ID', 'GT', 'chr') %in% colnames(markers)))
    stop('Missing cols in marker set')
  if (!all(c('sample_ID', 'pheno') %in% colnames(pheno)))
    stop('Missing cols in samples set')
  # NB, GT is 0 for normal, 1 for deletions, and 2 for duplications
  # NB, pheno is 1 for cases and 0 for controls

  # minimal table to work with
  dt <- merge(markers[, .(sample_ID, marker_ID, GT, chr)],
              pheno[, .(sample_ID, pheno)], by = 'sample_ID', all.x = T)
  # total number of cases and controls
  tot_cases <- pheno[pheno == 1, .N]
  tot_controls <- pheno[pheno == 0, .N]

  out <- data.table()
  # Run per chromosome
  for (cc in markers[, unique(chr)]) {
    message('Chromosome ', cc)
    dt_tmp <- unique(dt[chr == cc, ])
    dt_tmp_unique <- unique(dt_tmp[, .(marker_ID, GT)])
    markers_cc <- dt_tmp_unique[, marker_ID]
    gt_cc <- dt_tmp_unique[, GT]

    # fisher test
    if (test == 'fisher') {
      # function for mapply
      f_test_m <- function(marker, gt) {
        tmp <- dt_tmp[marker_ID == marker & GT == gt, ]

        if (tmp[, .N] > 0) {
          n_cases <- tmp[pheno == 1, .N]
          n_controls <- tmp[pheno == 0, .N]
          ft_matrix <- matrix(c(n_cases, tot_cases - n_cases,
                                n_controls, tot_controls - n_controls), ncol = 2)
          return(signif(as.numeric(fisher.test(ft_matrix)$p.value), digits = 2))
        } else return(NA)
      }
      # mapply
      ft_pval <- mapply(f_test_m, marker = markers_cc, gt = gt_cc)
      dt_tmp_unique[, pval := ft_pval]
    }
    out <- rbind(out, dt_tmp_unique)
  }
  return(out)
}


exclude_fixed_loci <- function(assoc_res, loci, boundary = 200000) {

  dt_f <- data.table()
  for (i in 1:loci[,.N]) {
    loc <- loci[i]
    dt_f <- rbind(dt_f, assoc_res[chr == loc$chr & start <= loc$end & end >= loc$start, ])
  }
  dt_no_f <- fsetdiff(assoc_res, dt_f)
  return(list(dt_f, dt_no_f))
}
