# README


## How to install

You can install the package via:

```
devtools::install_github("sinomem/CNValidatron_fl")
```

The package depends on two  non standard R packages,
R torch (https://torch.mlverse.org/) and BiocParallel
(https://www.bioconductor.org/packages/release/bioc/html/BiocParallel.html).
You might want to install them manually beforehand.


## Citation

If you use this software please cite the following publications:
https://doi.org/10.1101/2024.09.09.612035


## How to run

### Requirements

To perform the validation of a set of CNVs you will need the following files:

- tabix indexed intensity file for each sample
- SNPs file (preferably filtered)
- samples file (linking each sample to an intensity file)
- CNVs table

In general you will need all the files described in the our CNV calling protocol
(Montalbano et al., 2022, Current Protocol,
https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.621).
If you are interested in genome wide CNVs rather than in specific
loci (as in this case), you can just skip that section of the protocol.

Some example data in the correct formats is available in `data`.


### Run the prediction model

With that, you can run the actual model on your CNV table.
The main step of the pipeline run per sample. Here I'm showing a possible how to
run the program in batches but entirely on R. The for loop can easily be parallelised
outside R (e.g. across multiple jobs on a cluster).

```
# Pipeline #

library(CNValidatron)
setwd('path/to/working/directory')

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

  save_pngs_prediction(pngs_pt_batch, batch_cnvs, batch_samps, snps, no_parall = F)

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
```

Assuming the formats were correct and all went well, the resulting `all_preds`
data.table will contain the additional column `prob`. This is the prediction
probability of the CNV being true (between 0 and 1) for the specific GT, 1 for
deletions and 2 for duplications. The parameter `clean_out` can be set to `F` in
`make_predictions()` if you want to keep all three prediction probabilities
(`prob_true_del`, `prob_true_dup`, `prob_false`).

**NB**: `sample_ID` is used in the PNGs file name,
thus is best to not have special character in it (especially '/').


### Suggestions

If you are applying the program to a new dataset for the first time
I would strongly recommend you do some accuracy or at least precision
tests. If you don't have an in house solution to perform visual
validation you can use mine, https://github.com/SinomeM/shinyCNV.
It requires the same files and format as all the other CNV
packages I created so you should already have everything.


## How to train your own model

You might want to train a custom model for multiple reasons.
Your dataset might be very different for the ones I have access
to or you use case might be different (e.g. you are interested in
chromosomal abnormalities).


### Base case

In the base case (CNVs but in a very different dataset),
train a new custom model is fairly easy
assuming you already have a relatively big set of validated examples.

**NB**: This assumes you have a set of validated CNVs in
the following three classes: false, true deletion, true duplication.
The format should be the same as the described here
https://github.com/SinomeM/shinyCNV. That is column `vo`
for the human evaluation results (1: true, 2: false, 3: unknown)
and column `GT` for genotype (1: deletions, 2: duplications).   
Unknown CNVs will be ignored by the program.

```
# set BiocParall workers (for PNG generation)
options(MulticoreParam = MulticoreParam(workers = 16))

# load data
cnvs <- fread('path/to/validated_cnvs.txt')
samples <- fread('path/to/samples.txt')
snps <- fread('/path/to/snppos_filtered.txt')

# Split train and valitation 75/25
train_dt <- cnvs %>% group_by(GT) %>% sample_frac(size = 0.75)
train_dt <- as.data.table(train_dt)
valid_dt <- fsetdiff(cnvs, train_dt)
train_pt <- '/path/to/train'
valid_pt <- '/path/to/valid'

# save pngs
save_pngs_dataset(train_pt, train_dt, samples, snps, flip_chance = 0.5)
save_pngs_dataset(valid_pt, valid_dt, samples, snps, flip_chance = 0.5)

# training and validation dataloaders
bs <- 64
train_dl <- dataloader(torchvision::image_folder_dataset(
                         train_pt, transform = . %>% transform_to_tensor()),
                       batch_size = bs, shuffle = TRUE)
valid_dl <- dataloader(torchvision::image_folder_dataset(
                         valid_pt, transform = . %>% transform_to_tensor()),
                       batch_size = bs, shuffle = TRUE)

# Training

## Define the model etc
model <- convnet_dropout_5_10 %>%
  setup(
    loss = nn_cross_entropy_loss(),
    optimizer = optim_adam,
    metrics = list(
      luz_metric_accuracy()
    )
  )

# Learning rate finder
rates_and_losses_do <- model %>% lr_finder(train_dl, end_lr = 0.3)
rates_and_losses_do %>% plot() # REMEMBER TO CHECK AND UPDATE max_lr !!!!!

# Fit
fitted_dropout_5_10 <- model %>%
  fit(train_dl, epochs = 50, valid_data = valid_dl,
      callbacks = list(
        luz_callback_early_stopping(patience = 2),
        luz_callback_lr_scheduler(
          lr_one_cycle,
          max_lr = 0.01,
          epochs = 50,
          steps_per_epoch = length(train_dl),
          call_on = "on_batch_end"),
        luz_callback_model_checkpoint(path = "dropout_5_10_ukb_decode/"),
        luz_callback_csv_logger("dropout_5_10_ukb_decode.csv")
        ),
      verbose = TRUE)
# save the trained model
luz_save(fitted_dropout_5_10,
         'path/to/dropout_5_10_ukb_decode.rds')
```


### Other cases

If you want to use a different number of classes you need to change
at least the model definition in `R/.cnn_model.R`.

If you use case needs a different image, more/less pixels, different
set of tracks, different zoom etc, you might need to fork the repo
and change the function `plot_cnv()`.

In general for any complex use case, feel free to contact me
if that's the case and I'll see what I can do.

Unfortunately at the moment fine tuning of a pre-trained model
is not available, if you are an expert on R-torch and know how to
implement it, let me know!


## Pre-trained model

The model described in the manuscript is available here
https://doi.org/10.5281/zenodo.17174637.


## Bugs, Feature request, Collaborations and Contributions

Feel free to open an issue here on GitHub.
For collaboration and big extensions of the program, contact the
corresponding author of one of the mentioned papers.

Keep in mind I'm not a programmer and my main interest is
applying methods to do new research so be patient ;)

