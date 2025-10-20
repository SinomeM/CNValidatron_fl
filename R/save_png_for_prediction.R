#' Save the PNG Plots of CNVs for Prediction
#'
#' This function is used to create the dataset of PNG images for prediction
#'
#' @param root root folder for the dataset. Must not exists.
#' @param cnvs cnv data.table in the usual format
#' @param samps sample list in usual format
#' @param snps snps in the usual format
#' @param shrink_lrr see load_snps_tbx(), should be the same the model used in training
#'
#' @export
#'
#' @import data.table

save_pngs_prediction <- function(root, cnvs, samps, snps, shrink_lrr = 0.2,
                                 no_parall = F) {

  if (dir.exists(root)) warning('Root folder already exists!')
  dir.create(root, showWarning = F)

  # Create folders for each sample
  for (samp in cnvs[, unique(sample_ID)]) {
    samp_dir <- paste0(root, '/', samp, '/new/')
    dir.create(samp_dir, showWarnings = F, recursive = T)
  }

  FUN <- function(samp_id) {
    # Get all CNVs for this sample
    sample_cnvs <- cnvs[sample_ID == samp_id]
    sample_info <- samps[sample_ID == samp_id]
    
    # Read all tabix data for chromosomes needed by this sample
    chroms_needed <- sample_cnvs[, unique(chr)]
    tabix_data <- list()
    
    for (chr in chroms_needed) {
      tabix_data[[as.character(chr)]] <- load_snps_tbx(data.table(chr = chr),
                                                       sample_info, snps)
    }
    
    # Process each CNV for this sample
    for (i in 1:nrow(sample_cnvs)) {
      a <- sample_cnvs[i]
      
      # Get pre-loaded tabix data for this chromosome
      chr_data <- tabix_data[[as.character(a$chr)]]
      
      dt <- plot_cnv(a, sample_info, shrink_lrr = shrink_lrr, tabix_data = chr_data)
      n_real_snps <- dt[[2]]
      dt <- dt[[1]]
      
      pt <- paste0(root, '/', a$sample_ID, '/new/chr', a$chr, 
                   '_st', a$start, '_nsnp', n_real_snps, '.png')
      
      tryCatch({
        imager::save.image(imager::as.cimg(dt), pt)
      }, error = function(e) {
        print(paste(pt, 'failed.'))
      })
    }
    
    gc()
  }

  sample_ids <- cnvs[, unique(sample_ID)]
  
  if (no_parall)  null <- lapply(sample_ids, FUN)
  else null <- BiocParallel::bplapply(sample_ids, FUN)
}
