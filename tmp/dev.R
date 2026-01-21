devtools::build_vignettes()
devtools::document()
devtools::test()

devtools::check()
devtools::build()
devtools::install()

devtools::load_all()


# Minimal test run #

devtools::load_all()

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
preds[]

# One single false CNV
make_predictions(luz::luz_load('./joint.rds'),
                 pngs_pt, cnvs[1], return_pred_dt = F)[]

# One single true CNV
make_predictions(luz::luz_load('./joint.rds'),
                 pngs_pt, cnvs[3], return_pred_dt = F)[]