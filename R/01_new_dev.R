## (RE)START FROM SCRATCH THE DESIGN ##

setwd('~/Documents/CNValidatron_fl')
library(data.table)
library(ggplot2)

# Load CNV and samples table from UKB export
cnvs <- fread('~/Documents/UKB_data/edited/cvns.txt')
samples <- fread('~/Documents/UKB_data/edited/samples.txt')
loci <- fread('~/Documents/UKB_data/edited/loci.txt')
snps <- fread('~/Documents/UKB_data/snppos_filtered.txt')
cnvs <- cnvs[numsnp >= 20,]

# load visual validations
tmp <- list.files('~/Documents/UKB_data/visual_inspection2/')
tmp <- tmp[grep('visual', tmp)]
vi_cnv <- data.table()
for (i in tmp)
  vi_cnv <- rbind(vi_cnv,
                  fread(paste0('~/Documents/UKB_data/visual_inspection2/', i)))

# check for duplicates
vi_cnv[, ':=' (V1 = NULL, index = NULL)]
nrow(unique(vi_cnv)) == nrow(vi_cnv)
nrow(unique(vi_cnv[, .(sample_ID, chr, start, GT)])) == nrow(vi_cnv)
tmp <- vi_cnv[duplicated(vi_cnv[, .(sample_ID, chr, start, GT)])]
vi_cnv[sample_ID %in% tmp$sample_ID & chr %in% tmp$chr & start %in% tmp$start &
       GT %in% tmp$GT, ]

# remember to deal with the modified boundaries
vi_cnv[Visual_Output == -8, ]

# It is useful that here each locus is unique
length(unique(cnvs$locus)) == nrow(cnvs)
# --------------------------------------------------------------------------- #




## STEP 1 ##
# produce the plots in a format the algorithm can undrstand best

# option 1: pixelated image of classic BAF LRR in a square figure
ii <- 12
cc <- vi_cnv[numsnp > 40 & Visual_Output == 1, ][ii]

ior <- 3; shlrr <- 0.25; width <- 50

print(check_cnv(cc, samples[sample_ID == cc[, sample_ID], ], snps = snps,
                in_out_ratio = ior, shrink_lrr = shlrr))

tmp <- plot_cnv(cc, samples[sample_ID == cc[, sample_ID], ], snps = snps,
                tmp_plot = 2, w = width, in_out_ratio = ior, shrink_lrr = shlrr)
# --------------------------------------------------------------------------- #



## STEP 2 ##
# create and store the data for each CNV

# - In torch the plots will be tensors
# - In torch the data handling is done via a dataloader()
# - Since the function torchvision::image_folder_dataset() exists, the
#   easiest thing to do is save all plots as PNG in the expected file
#   structure

library(torch)
library(torchvision)
library(torchdatasets)
library(luz)
source('./R/02_new_dev_functions.R')

# just to test data structure etc really, way too small to even early testing
# of an actual model
tmp <- vi_cnv[numsnp > 40 & Visual_Output %in% 1:3, ]
train_test <- tmp[sample(1:nrow(tmp), round(nrow(tmp)*0.7) )]
valid_test <- fsetdiff(tmp, train_test)

# run if necessary
train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/train'
if (F)
  save_pngs_dataset(train_pt, train_test, samples, snps)
valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/valid'
if (F)
  save_pngs_dataset(valid_pt, valid_test, samples, snps)
# by default the images are 64x64

tdt <- torchvision::image_folder_dataset(
  train_pt, transform = . %>% transform_to_tensor())
vdt <- torchvision::image_folder_dataset(
  valid_pt, transform = . %>% transform_to_tensor())

# Is this the correct format? Test on a simple model
train_dl <- dataloader(tdt, batch_size = 128, shuffle = TRUE)
valid_dl <- dataloader(vdt, batch_size = 128, shuffle = TRUE)

batch <- train_dl %>%
  dataloader_make_iter() %>%
  dataloader_next()

dim(batch$x)
batch$y

convnet <- nn_module(
  "convnet",
  initialize = function() {
    self$features <- nn_sequential(
      nn_conv2d(3, 64, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(64, 128, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(128, 256, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(256, 512, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(512, 1024, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_adaptive_avg_pool2d(c(1, 1))
    )
    self$classifier <- nn_sequential(
      nn_linear(1024, 1024),
      nn_relu(),
      nn_linear(1024, 1024),
      nn_relu(),
      nn_linear(1024, 5)
    )
  },
  forward = function(x) {
    x <- self$features(x)$squeeze()
    x <- self$classifier(x)
    x
  }
)

fitted <- convnet %>%
  setup(
    loss = nn_cross_entropy_loss(),
    optimizer = optim_adam,
    metrics = list(
      luz_metric_accuracy()
    )
  ) %>%
  fit(train_dl,
      epochs = 50,
      valid_data = valid_dl,
      verbose = TRUE
      )
# --------------------------------------------------------------------------- #



## STEP 3 ##
# data augmentation

# - torch includes functions to perform data augmentation quite easily
#   but since the task is quite weird I think I will do at least some
#   manually. E.g. a lot of plots have big holes in them or one side
#   almost completely missing, this is something I can play with.
# - Image flipping on the center vertical axis is something I can leave
#   to torch since it's trivial, rotations does not make any sense.
# - Other? one example could be "transplant" the CNV to a different region,
#   not sure how much would make sense

if (F) {
  tdt <- torchvision::image_folder_dataset(
    train_pt, transform = . %>% transform_to_tensor() %>%
                                transform_random_horizontal_flip())
  vdt <- torchvision::image_folder_dataset(
    valid_pt, transform = . %>% transform_to_tensor() %>%
                                transform_random_horizontal_flip())
}
# I'm not sure if transform_random_horizontal_flip() adds a flipped copy
# or flip the actual example. In this context I think I want the first, since
# each example is quite expensive to produce. Since it's not clear in the
# documentation, I think I'll do it manually.

# it is implemented in save_pngs_dataset()
source('./R/02_new_dev_functions.R')

# run if necessary
if (F) {
  tmp <- vi_cnv[numsnp > 40 & Visual_Output %in% 1:3, ]
  train_test <- tmp[sample(1:nrow(tmp), round(nrow(tmp)*0.7) )]
  valid_test <- fsetdiff(tmp, train_test)

  npx = 64
  train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/train'
  save_pngs_dataset(train_pt, train_test, samples, snps, w = npx, flip_chance = 0.7)
  valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/valid'
  save_pngs_dataset(valid_pt, valid_test, samples, snps, w = npx, flip_chance = 0.7)
}

# the second layer (regarding the 'holes' in the data) is a bit more tricky
# and I'll work on it in the future

# --- --- --- ---

# Finally "transplanting" the CNV in a different region (possibly even
# a different sample) is something I'm not sure it's helpful, I'll consider it
# in the future

# --- --- --- ----


# --------------------------------------------------------------------------- #



## STEP 4 ##
# Is the data in the best shape? Can it improve?

# At the moment I think it is worth trying as it is, I need more examples
# and a proper model, anyway I'll think again about it in the future.

# --- --- --- ---

# --------------------------------------------------------------------------- #



## STEP 5 ##
# Build a proper model not just for early testing
source('./R/02_new_dev_functions.R')

train_dl <- dataloader(torchvision::image_folder_dataset(
                         train_pt, transform = . %>% transform_to_tensor()),
                       batch_size = 128, shuffle = TRUE)
valid_dl <- dataloader(torchvision::image_folder_dataset(
                         valid_pt, transform = . %>% transform_to_tensor()),
                       batch_size = 128, shuffle = TRUE)

npx
classes <- 5 # T del/dup, U del/dup, F

convnet_dropout <- nn_module(
  "convnet",
  initialize = function() {
    self$features <- nn_sequential(
      nn_conv2d(3, npx, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx, npx*2, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*2, npx*4, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*4, npx*8, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*8, npx*16, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_adaptive_avg_pool2d(c(1, 1)),
      nn_dropout2d(p = 0.05),
    )
    self$classifier <- nn_sequential(
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, classes)
    )
  },
  forward = function(x) {
    x <- self$features(x)$squeeze()
    x <- self$classifier(x)
    x
  }
)

convnet_batchnorm <- nn_module(
  "convnet",
  initialize = function() {
    self$features <- nn_sequential(
      nn_conv2d(3, npx, kernel_size = 3, padding = 1),
      nn_batch_norm2d(npx),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(npx, npx*2, kernel_size = 3, padding = 1),
      nn_batch_norm2d(npx*2),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(npx*2, npx*4, kernel_size = 3, padding = 1),
      nn_batch_norm2d(npx*4),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(npx*4, npx*8, kernel_size = 3, padding = 1),
      nn_batch_norm2d(npx*8),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_conv2d(npx*8, npx*16, kernel_size = 3, padding = 1),
      nn_batch_norm2d(npx*16),
      nn_relu(),
      nn_adaptive_avg_pool2d(c(1, 1)),
    )
    self$classifier <- nn_sequential(
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_batch_norm1d(npx*16),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_batch_norm1d(npx*16),
      nn_linear(npx*16, classes)
    )
  },
  forward = function(x) {
    x <- self$features(x)$squeeze()
    x <- self$classifier(x)
    x
  }
)

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
        luz_callback_early_stopping(patience = 5),
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

# same as the other one
model <- convnet_batchnorm %>%
  setup(
    loss = nn_cross_entropy_loss(),
    optimizer = optim_adam,
    metrics = list(
      luz_metric_accuracy()
    )
  )
rates_and_losses_bn <- model %>% lr_finder(train_dl, end_lr = 0.3)
rates_and_losses_bn %>% plot() # REMEMBER TO CHECK AND UPDATE max_lr !!!!!

fitted_batchnorm <- model %>%
  fit(train_dl, epochs = 50, valid_data = valid_dl,
      callbacks = list(
        luz_callback_early_stopping(patience = 5),
        luz_callback_lr_scheduler(
          lr_one_cycle,
          max_lr = 0.00001,
          epochs = 50,
          steps_per_epoch = length(train_dl),
          call_on = "on_batch_end"),
        luz_callback_model_checkpoint(path = "cpt_batchnorm/"),
        luz_callback_csv_logger("logs_batchnorm.csv")
        ),
      verbose = TRUE)
# even at low max_lr it seems this model tends to overfit at the moment,
# maybe it will be more helpful when I have a larger dataset


# --------------------------------------------------------------------------- #




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
  train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/test3/train'
  valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/test3/valid'
  save_pngs_dataset(train_pt, train_test, samples, snps, w = npx, flip_chance = 0.7)
  save_pngs_dataset(valid_pt, valid_test, samples, snps, w = npx, flip_chance = 0.7)
}

if (F) {
  # minsnp 35 or lower
  tmp <- vi_cnv[numsnp > 35 & Visual_Output %in% 1:3, ]
  train_test <- tmp[sample(1:nrow(tmp), round(nrow(tmp)*0.75) )]
  valid_test <- fsetdiff(tmp, train_test)
  train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/test2/train'
  valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/test2/valid'
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

batch <- train_dl %>%
  dataloader_make_iter() %>%
  dataloader_next()

dim(batch$x)
batch$y

npx
classes <- 5 # T del/dup, U del/dup, F

# at the moment this seems the best easy model
convnet_dropout <- nn_module(
  "convnet",
  initialize = function() {
    self$features <- nn_sequential(
      nn_conv2d(3, npx, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx, npx*2, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*2, npx*4, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*4, npx*8, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_max_pool2d(kernel_size = 2),
      nn_dropout2d(p = 0.05),
      nn_conv2d(npx*8, npx*16, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_adaptive_avg_pool2d(c(1, 1)),
      nn_dropout2d(p = 0.05),
    )
    self$classifier <- nn_sequential(
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, npx*16),
      nn_relu(),
      nn_dropout(p = 0.05),
      nn_linear(npx*16, classes)
    )
  },
  forward = function(x) {
    x <- self$features(x)$squeeze()
    x <- self$classifier(x)
    x
  }
)

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
        luz_callback_early_stopping(patience = 5),
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
}

cnvs
tmp <- copy(vi_cnv)
tmp[, ':=' (V1 = NULL, index = NULL, Visual_Output = NULL)]

test_dt <- fsetdiff(cnvs, tmp)

if (F) {
  for (i in 1:nrow(test_dt)) {
    print(i)
    dt <- load_snps_tbx(test_dt[i], samples, snps)
    if (nrow(dt) < 40) test_dt[i, rm := T]
  }
  test_dt <- test_dt[is.na(rm), ]
  saveRDS(test_dt, 'tmp/test_dt.rds')
}

test_dt <- readRDS('tmp/test_dt.rds')
test_dt[, ':=' (Visual_Output = 5, rm = NULL)]


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
