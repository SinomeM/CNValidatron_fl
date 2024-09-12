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

- Regarding CNV data handling and generation:
  https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.621
- Regarding the actual model:
  insert_here


## How to run

To perform the validation of a set of CNVs you will need the following files:

- tabix indexed intensity file for each sample
- SNPs file (preferably filtered)
- samples file (linking each sample to an intensity file)
- CNVs table

In general you will need all the files described in the CNV calling protocol
mentioned above (Montalbano et al., 2022, Current Protocol).
If you are interested in genome wide CNVs rather than in specific
loci, you can just skip that section of the protocol.

With that, you can run the actual model on your CNV table as shown in
the following code snippet:

```
# load necessary objects
snps <- fread('/path/to/snppos_filtered.txt')
cnvs <- fread('/path/toc/nvs_filtered.txt')
samples <- fread('/path/tol/samples_list.txt')

# select the folder for PNG files
png_pt <- '/path/to/folder'

# set BiocParall parallel worker limit
BiocParallel::register(BiocParallel::MulticoreParam(workers=2))

# save the PNGs for all CNVs
save_pngs_prediction(pred_pt, cnvs, samples, snps)

# run the prediction algoritm
preds <- make_predictions(luz::luz_load('/path/to/dropout_5_10_ukb_decode.rds'),
                          png_pt, cnvs)

# select predicted true CNVs with probability above 0.75
true_cnvs <- pred[pred %in% 2:3 & pred_prob >= 0.75, ]

# the model has three categories
# 1: False
# 2: True Deletion
# 3: True Duplication
```

Notice that at the moment the sample_ID is used in the PNGs file name,
thus is best to not have special character in it (especially '/').


## How to train your on model



## Pre-trained model

The model described in the paper above is available upon request.


## Bugs, Feature request, Collaborations and Contributions

Feel free to open an issue here on GitHub.
For collaboration and big extensions of the program, contact the 
corresponding author in one of the mentioned papers.
