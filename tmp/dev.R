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
save_pngs_prediction(pngs_pt, cnvs, samples, snps, no_parall = F)

# Run the prediction algoritm
preds <- make_predictions(luz::luz_load('./joint.rds'),
                          pngs_pt, cnvs, return_pred_dt = F)
preds[]

# One single false CNV
make_predictions(luz::luz_load('./joint.rds'),
                 pngs_pt, cnvs[1], return_pred_dt = F)[]

# One single true CNV
make_predictions(luz::luz_load('./joint.rds'),
                 pngs_pt, cnvs[3], return_pred_dt = F)[]




# Test one single CNV #
devtools::load_all()
test_cnv <- cnvs[3]
unlink(pngs_pt, recursive = TRUE)
save_pngs_prediction(pngs_pt, test_cnv, samples, snps, no_parall = F)

# Run the prediction algoritm
test <- make_predictions(luz::luz_load('./joint.rds'),
                          pngs_pt, test_cnv, return_pred_dt = T)
test[] # good
