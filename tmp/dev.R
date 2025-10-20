devtools::build_vignettes()
devtools::document()
devtools::test()

devtools::check()
devtools::build()
devtools::install()

devtools::load_all()


# Minimal test run #

# load necessary objects
snps <- fread('./data/hd_1kG_hg19.snppos.filtered.test.gz')
cnvs <- fread('./data/cnvs.txt')
cnvs[, prob := NULL]
samples <- fread('./data/samples_list.txt')

# select the folder for PNG files
pngs_pt <- './tmp/pngs'

# set BiocParall parallel worker limit
BiocParallel::register(BiocParallel::MulticoreParam(workers=2))

# Save the PNGs for all CNVs
unlink(pngs_pt, recursive = TRUE)
save_pngs_prediction(pngs_pt, cnvs[chr != 22, ], samples, snps, no_parall = F)
# data for chromosome 22 seems to be missing from the tabix (???)

# Run the prediction algoritm
devtools::load_all()
preds <- make_predictions(luz::luz_load('./joint.rds'),
                          pngs_pt, cnvs, return_pred_dt = F)
preds


# Pipeline #

devtools::load_all()
# data 
snps <- fread('./data/hd_1kG_hg19.snppos.filtered.test.gz')
cnvs <- fread('./data/cnvs.txt')
cnvs <- cnvs[, prob := NULL][chr != 22, ]
samples <- fread('./data/samples_list.txt')

# select the folder for PNG files
pngs_pt <- './tmp/pngs'

# set BiocParall parallel worker limit
BiocParallel::register(BiocParallel::MulticoreParam(workers=2))

# Create batches, these can be one per sample or multiple samples per batch.
batches <- 2
samples[, batch := sample(1:batches, .N, replace = T)]

# Batches can be parallelized outside of R if needed.
for (b in 1:batches) {
  # select samples and CNVs for this batch
  batch_samps <- samples[batch == b, sample_ID]
  batch_cnvs <- cnvs[sample_ID %in% batch_samps, ]

  # batch subfolder
  pngs_pt_batch <- paste0(pngs_pt, '/batch_', b)
  # clean just to be sure
  unlink(pngs_pt_batch, recursive = TRUE)

  save_pngs_prediction(pngs_pt_batch, batch_cnvs, samples, snps, no_parall = F)

  preds <- make_predictions(luz::luz_load('./joint.rds'),
                            pngs_pt_batch, batch_cnvs, return_pred_dt = F)
  # remove PNGs after prediction has completed
  unlink(pngs_pt_batch, recursive = TRUE)

  # save predictions for this batch
  dir.create('./tmp', showWarnings = F) # create tmp folder if needed
  fwrite(preds, paste0('./tmp/preds_batch_', b, '.txt'), sep = '\t')
}

# remove main PNG folder
unlink(pngs_pt, recursive = TRUE)

# Combine all batch predictions
all_preds <- data.table()
for (b in 1:batches) {
  batch_preds <- fread(paste0('./tmp/preds_batch_', b, '.txt'))
  all_preds <- rbind(all_preds, batch_preds)
  # remove batch prediction files
  file.remove(paste0('./tmp/preds_batch_', b, '.txt'))
}

all_preds
