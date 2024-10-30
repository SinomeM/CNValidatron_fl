#' Quick fisher test association function
#'
#' Simple function to perform a fisher test for each
#' marker against a phenotype.
#'
#'
#'
#'
#' @import data.table
#'
#' @export


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
  return(list(in_loci = dt_f, not_in_loci = dt_no_f))
}



#' Process results from quick assoc
#'
#' Process `quick_assoc()` outputs for easier visualization.
#' Perform p value adjustment, filter based on max p value
#' and exclude fixed loci.
#'
#'
#' @import data.table
#'
#' @export

process_quick_res <- function(quick_res, markers, padj = F, max_pval = 0.1,
                              loci, loci_boundaries = 250000) {
  if (padj) {
    quick_res[, adj_pval := signif(p.adjust(pval, method = 'fdr'), 3)]
    quick_res <- quick_res[adj_pval <= max_pval,]
  }
  else quick_res <- quick_res[pval <= max_pval,]

  dt <- merge(quick_res, unique(markers[, .(marker_ID, chr, start, end)]),
              by = 'marker_ID')
  res <- CNValidatron:::exclude_fixed_loci(dt, loci, 250000)
}




#' Run logisti regression on a set of CNV markers
#'
#'
#'
#'
#' @import data.table
#'
#' @export

run_assoc_logistic <- function(markers, scan_res, pheno, padj = F,
                               max_pval = 0.1, min_carriers = 5) {

  # re unique to be always sure
  markers <- unique(markers)
  scan_res <- unique(scan_res)
  pheno <- unique(pheno)

  # select markers
  if (padj)
    dt <- markers[marker_ID %in% scan_res[adj_pval <= max_pval, marker_ID], ]
  else
    dt <- markers[marker_ID %in% scan_res[pval <= max_pval, marker_ID], ]

  all_markers <- dt[, unique(marker_ID)]
  l_m <- length(all_markers)

  # run logistic regression for each selected marker
  res <- data.table()
  for (i in 1:l_m) {
    message(i, ' of ', l_m, ' markers')
    mm <- all_markers[i]
    dtm <- dt[marker_ID == mm, ]

    # number of carriers
    n_del <- dtm[GT == 1, .N]
    n_dup <- dtm[GT == 2, .N]
    # if none above 5 skip marker
    if (n_del < min_carriers & n_dup < min_carriers) next

    # pheno file for the marker
    tmp <- pheno[, .(sample_ID, pheno, age, gender)]
    if (n_del < min_carriers)
      tmp[sample_ID %in% dtm[GT == 1, sample_ID], var := 1]
    if (n_dup < min_carriers)
      tmp[sample_ID %in% dtm[GT == 2, sample_ID], var := 2]
    tmp[is.na(var), var := 0]

    # HERE DUP carriers will be treated as controls if the DEL has
    # more than 5 but the DUP not. In GDK is done differently

    # run association
    mod <- gam(pheno ~ var + s(age) + gender, data = tmp, family = binomial)

    # format results
    pt <- summary(mod)$p.table
    tmp <- cbind(rownames(pt), as.data.table(pt))
    colnames(tmp) <- c('var', 'log_OR', 'SE', 'Z', 'pval')
    if (n_del >= 5)
      res <- rbind(res, data.table(marker_ID = mm, GT = 1, n<-car = n_del,
                                   log_OR = signif(tmp[var == 'var1', log_OR], 3),
                                   SE = signif(tmp[var == 'var1', SE], 3),
                                   pval = signif(tmp[var == 'var1', pval], 3)))
    if (n_dup >= 5)
      res <- rbind(res, data.table(marker_ID = mm, GT = 2, n<-car = n_del,
                                   log_OR = signif(tmp[var == 'var2', log_OR], 3),
                                   SE = signif(tmp[var == 'var2', SE], 3),
                                   pval = signif(tmp[var == 'var2', pval], 3)))

  }
  return(list(res, markers, pheno))
}
