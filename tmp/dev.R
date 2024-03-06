
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


library(data.table)
samps <- fread('~/Documents/CNValidatron_trained_models/second_iteration/all_samples_15k.txt')
cnvs <- fread('~/Documents/CNValidatron_trained_models/second_iteration/all_training_example_8k.txt')
a <- cnvs[1]
b <- samps[sample_ID == a$sample_ID, ]
plot_cnv(a, b)



library(data.table)
devtools::load_all()

cnvs <- fread('../../UKB_GW_CNVs/cnvs_pred.txt')[pred %in% 2:3 & pred_prob >= 0.9, ]
# mild filters on outliers
cnvs <- cnvs[length <= 10000000, ]
cnvs <- cnvs[!sample_ID %in% cnvs[ , .N, by = sample_ID][N >10, sample_ID], ]
fwrite(cnvs[, .(chr, start, end, GT)], '../../UKB_GW_CNVs/igv/cnvs.bed', col.names = F, sep = '\t')

