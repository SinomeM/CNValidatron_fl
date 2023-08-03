
## Model training start to end ##

# All lines to test the model(s) training should be here
if (T) {
  setwd('~/Documents/CNValidatron_fl')
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
  library(torch)
  library(torchvision)
  library(luz)
  source('./R/02_new_dev_functions.R')
  options(MulticoreParam = MulticoreParam(workers = 12))

  # Load CNV and samples table from UKB export
  samples <- fread('~/Documents/UKB_data/edited/samples.txt')
  snps <- fread('~/Documents/UKB_data/snppos_filtered.txt')
  # load visual validations
  tmp <- list.files('~/Documents/UKB_data/visual_inspection2/')
  tmp <- tmp[grep('visual', tmp)]
  vi_cnv <- data.table()
  for (i in tmp)
    vi_cnv <- rbind(
      vi_cnv,fread(paste0('~/Documents/UKB_data/visual_inspection2/', i)))
}

# only larger CNVs for now
tmp <- vi_cnv[numsnp > 40 & Visual_Output %in% 1:3, ]
train_test <- tmp[sample(1:nrow(tmp), round(nrow(tmp)*0.80) )]
valid_test <- fsetdiff(tmp, train_test)

# class imbalance?
train_test[, .N, by = c('Visual_Output', 'GT')]
valid_test[, .N, by = c('Visual_Output', 'GT')]
# there is some class imbalance, however this is due to the true distribution
# of CNV calls produced by PennCNV so in a way it should be this way (?)

# save plots
npx = 64
if (F) {
  # minsnp 40, default
  train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/min40snps/train'
  valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/min40snps/valid'
  save_pngs_dataset(train_pt, train_test, samples, snps, w = npx, flip_chance = 0.7)
  save_pngs_dataset(valid_pt, valid_test, samples, snps, w = npx, flip_chance = 0.7)
}

if (F) {
  # minsnp 35 or lower
  tmp <- vi_cnv[numsnp > 35 & Visual_Output %in% 1:3, ]
  train_test <- tmp[sample(1:nrow(tmp), round(nrow(tmp)*0.75) )]
  valid_test <- fsetdiff(tmp, train_test)
  train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/min35snps/train'
  valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/min35snps/valid'
  save_pngs_dataset(train_pt, train_test, samples, snps, w = npx, flip_chance = 0.7)
  save_pngs_dataset(valid_pt, valid_test, samples, snps, w = npx, flip_chance = 0.7)
}

# training and validation dataloaders
bs <- 64
train_dl <- dataloader(torchvision::image_folder_dataset(
                         train_pt, transform = . %>% transform_to_tensor()),
                       batch_size = bs, shuffle = TRUE)
valid_dl <- dataloader(torchvision::image_folder_dataset(
                         valid_pt, transform = . %>% transform_to_tensor()),
                       batch_size = bs, shuffle = TRUE)

npx
classes <- 5 # T del/dup, U del/dup, F

# learning rate finder
model <- convnet_dropout %>%
  setup(
    loss = nn_cross_entropy_loss(),
    optimizer = optim_adam,
    metrics = list(
      luz_metric_accuracy()
    )
  )
rates_and_losses_do <- model %>% lr_finder(train_dl, end_lr = 0.3)
rates_and_losses_do %>% plot() # REMEMBER TO CHECK AND UPDATE max_lr !!!!!

# fit
fitted_dropout <- model %>%
  fit(train_dl, epochs = 50, valid_data = valid_dl,
      callbacks = list(
        luz_callback_early_stopping(patience = 2),
        luz_callback_lr_scheduler(
          lr_one_cycle,
          max_lr = 0.01,
          epochs = 50,
          steps_per_epoch = length(train_dl),
          call_on = "on_batch_end"),
        luz_callback_model_checkpoint(path = "cpt_dropout/"),
        luz_callback_csv_logger("logs_dropout.csv")
        ),
      verbose = TRUE)

# save the preliminary fitted model
luz_save(fitted_dropout, 'tmp/fitted_dropout.rds')

# --------------------------------------------------------------------------- #



## Test predictions and speedup new examples generation ##

# the models seems decent enough to be used in some testing. My plan
# is to use it to group the CNVs I have to validate so that it would
# be much easier to do (visual inspection is much quicker if the calls
# are similar).

# Tho models works best with CNVs of at least 40 SNPs so that will be my focus
# for the moment

if (T) {
  setwd('~/Documents/CNValidatron_fl')
  source('./R/02_new_dev_functions.R')
  library(BiocParallel)
  library(data.table)
  library(torch)
  library(torchvision)
  library(luz)
  fitted_model <- luz_load('tmp/fitted_dropout.rds')
  options(MulticoreParam = MulticoreParam(workers = 12))

  cnvs <- fread('~/Documents/UKB_data/edited/cvns.txt')
  samples <- fread('~/Documents/UKB_data/edited/samples.txt')
  snps <- fread('~/Documents/UKB_data/snppos_filtered.txt')
  cnvs <- cnvs[numsnp >= 40,]
  # load visual validations
  tmp <- list.files('~/Documents/UKB_data/visual_inspection2/')
  tmp <- tmp[grep('visual', tmp)]
  vi_cnv <- data.table()
  for (i in tmp)
    vi_cnv <- rbind(
      vi_cnv,fread(paste0('~/Documents/UKB_data/visual_inspection2/', i)))

  cnvs
  tmp <- copy(vi_cnv)
  tmp[, ':=' (V1 = NULL, index = NULL, Visual_Output = NULL)]
  test_dt <- fsetdiff(cnvs, tmp)
}

if (F) {
  for (i in 1:nrow(test_dt)) {
    print(i)
    dt <- load_snps_tbx(test_dt[i], samples, snps)
    if (nrow(dt) < 40) test_dt[i, rm := T]
  }
  test_dt <- test_dt[is.na(rm), ]

  for (i in 1:nrow(test_dt)) {
    print(i)
    a <- test_dt[i]
    pt <- paste0('tmp/unvalidated/', a$sample_ID, '_', a$start, '.png')

    dt <- plot_cnv(a, samples[sample_ID == a[, sample_ID], ], snps = snps,
                   w = 64, in_out_ratio = 3, shrink_lrr = 0.2)
    if (nrow(dt) == 0) {
      test_dt[i, rm := T]
      next
    }

    dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis
    imager::save.image(imager::as.cimg(dt), pt)
  }
  test_dt <- test_dt[is.na(rm), ]
  test_dt[, rm := NULL]
  saveRDS(test_dt, 'tmp/test_dt.rds')
}

test_dt <- readRDS('tmp/test_dt.rds')
uneval_dataset <- image_folder_dataset(
  'tmp/unvalidated/', transform = . %>% transform_to_tensor())
fitted_model
# the output is raw logits
uneval_pred <- predict(fitted_model, uneval_dataset)
# convert to probabilities
uneval_pred <- 1 / (1 + exp(-uneval_pred))


# the classes are
# false:   1
# tru_del: 2
# tru_dup: 3
# unk_del: 4
# unk_dup: 5


dt <- data.table(ix = uneval_dataset$samples[[1]],
                 class1 = round(as.numeric(uneval_pred[,1]), 5),
                 class2 = round(as.numeric(uneval_pred[,2]), 5),
                 class3 = round(as.numeric(uneval_pred[,3]), 5),
                 class4 = round(as.numeric(uneval_pred[,4]), 5),
                 class5 = round(as.numeric(uneval_pred[,5]), 5))

for (i in 1:nrow(dt))
  dt[i, pred := which.max(c(class1, class2, class3, class4, class5))]

dt[, sample_ID := gsub('/home/simone/Documents/CNValidatron_fl/tmp/unvalidated/tmp/', '', ix)]
dt[, start := gsub('\\d+_', '', sample_ID)][, start := as.integer(gsub('.png', '', start))]
dt[, sample_ID := as.integer(gsub('_\\w+.png', '', sample_ID))]
dt[, ix := NULL]


dt[, locus := paste0('locus', 1:nrow(dt))]
dt <- merge(dt, cnvs[, .(sample_ID, chr, start, end, numsnp,
                         length, type, conf, batch, GT, CN)],
            by = c('sample_ID', 'start'))
loci <- cnvs[, .(locus, chr, start, end)]

fwrite(loci, '../UKB_data/visual_inspection3/loci.txt', sep = '\t')
fwrite(loci[, .(locus)], '../UKB_data/visual_inspection3/loci_Sel.txt', col.names = F)
for (i in 1:5)
  fwrite(
    dt[pred == i, ][, .(sample_ID, chr, start, end, GT, CN, numsnp,
                        length, type, conf, locus, batch, pred)],
    paste0('../UKB_data/visual_inspection3/putative', i, '.txt'), sep = '\t')

