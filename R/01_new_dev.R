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
train_test <- vi_cnv[numsnp > 40 & Visual_Output %in% 1:3, ]
valid_test <- vi_cnv[] # TBD!!!!

# run if necessary
train_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/train'
if (F)
  save_pngs_dataset(train_pt, train_test, samples, snps)
valid_pt <- '/home/simone/Documents/CNValidatron_fl/tmp/valid'
if (F)
  save_pngs_dataset(valid_pt, valid_test, samples, snps)

tdt <- torchvision::image_folder_dataset(
  train_pt, transform = . %>% transform_to_tensor())
vdt <- torchvision::image_folder_dataset(
  valid_pt, transform = . %>% transform_to_tensor())

# Is this the correct format? Test on a simple model
train_dl <- dataloader(tdt, batch_size = 100, shuffle = TRUE)
valid_dl <- dataloader(vdt, batch_size = 100, shuffle = TRUE)
str(train_dl)

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
      nn_linear(1024, 200)
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
# data agumentation

# - torch includes functions to perform data agumentation quite easily
#   but since the task is quite wierd I think I will do at least some
#   manually. E.g. a lot of plots have big holes in them or one side
#   almost completely missing, this is something I can play with.
# - Image flipping on the center vertical axis is something I can leave
#   to torch since it's trivial, rotations does not make any sense.
# - Other? one example could be "transplant" the CNV to a different region,
#   not sure how much would make sense

# --------------------------------------------------------------------------- #



## STEP 4 ##
# Is the data in the best shape? Can it improve?

# --------------------------------------------------------------------------- #



## STEP 5 ##
# Build a proper model not just for early testing

# --------------------------------------------------------------------------- #
