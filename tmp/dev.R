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
debugonce(save_pngs_prediction)
save_pngs_prediction(pngs_pt, cnvs[sample_ID == 'sample2',], samples, snps, batches = 2, no_parall = T)
traceback()

# run the prediction algoritm
model_pt <- paste0('~/Desktop/pdrive/simone/01.Work/projects/',
                   'CNValidatron_models_2025/trained_models/joint.rds')
preds <- make_predictions(luz::luz_load(model_pt),
                          pngs_pt, cnvs)