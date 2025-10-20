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
samples <- fread('./data/samples_list.txt')

# select the folder for PNG files
pngs_pt <- './tmp/pngs'

# set BiocParall parallel worker limit
BiocParallel::register(BiocParallel::MulticoreParam(workers=2))

# save the PNGs for all CNVs
# debugonce(save_pngs_prediction)
devtools::load_all()
unlink(pngs_pt, recursive = TRUE)
save_pngs_prediction(pngs_pt, cnvs[chr != 22, ], samples, snps, no_parall = T)
traceback()
# data for chromosome 22 seems to be missing from the tabix (???)

# OK up to here


# run the prediction algoritm
preds <- make_predictions(luz::luz_load('./joint.rds'),
                          pngs_pt, cnvs)