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
cnvs[, prob := NULL]
samples <- fread('./data/samples_list.txt')

# select the folder for PNG files
pngs_pt <- './tmp/pngs'

# set BiocParall parallel worker limit
BiocParallel::register(BiocParallel::MulticoreParam(workers=2))

# Save the PNGs for all CNVs
unlink(pngs_pt, recursive = TRUE)
save_pngs_prediction(pngs_pt, cnvs[chr != 22, ], samples, snps, no_parall = F)
# data for chromosome 22 seems to be missing from the tabix (???)

# Run the prediction algoritm
devtools::load_all()
preds <- make_predictions(luz::luz_load('./joint.rds'),
                          pngs_pt, cnvs, return_pred_dt = F)
preds[]



# Grad-CAM visualization test #
devtools::load_all()

# Load your trained model
model <- luz::luz_load('./joint.rds')

# Select a sample PNG
img_path <- list.files('./tmp/pngs/sample1/new', pattern = '\\.png$',
                       full.names = TRUE)[7]

# Generate Grad-CAM
gradcam_result <- generate_gradcam(model, img_path,  target_class = 2)

# Plot the results
plot_gradcam(gradcam_result)

# Or save to file
plot_gradcam(gradcam_result, save_path = './tmp/gradcam_visualization.png')

# You can also specify a target class
gradcam_result_class0 <- generate_gradcam(model, img_path, target_class = 0)
plot_gradcam(gradcam_result_class0)
