
devtools::build_vignettes()

devtools::document()

devtools::test()

devtools::load_all()

devtools::check()

devtools::build()

devtools::install()

library(data.table)
pt <- '~/Documents/CNValidatron_fl/tmp/prova'
samps <- fread('~/Documents/UKB_data/edited/samples.txt')[sample(200),]
cnvs <- fread('~/Documents/CNValidatron_trained_models/all_visual_inspection_min40snps.txt')[
          sample_ID %in% samps$sample_ID,]
snps <- fread('~/Documents/UKB_data/snppos_filtered.txt')

save_pngs_prediction(pt, cnvs, samps, snps,  w = 64)


devtools::load_all()
pred <- make_predictions(luz::luz_load('~/Documents/CNValidatron_trained_models/fitted_dropout_all.rds'),
                         pt, cnvs)
pred

