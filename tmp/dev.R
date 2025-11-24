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



# Grad CAM plots on some examples #

# Find very good CNVs to use as examples, both deletions and duplications
library(data.table)
files <- gsub('\\.tabix.+', '', list.files('./data/1000Genomes_data/tbx/'))
samples <- fread('./data/1000Genomes_data/common_samples.txt')[hd_id %in% files, ]
samples <- samples[, .(sample_ID = Sample_Group,
            file_path_tabix = paste0('/Users/sm/Documents/GitHub/CNValidatron_fl/data/1000Genomes_data/tbx/',
            hd_id, '.tabix.gz'))]
dt <- fread('./data/1000Genomes_data/1kG_cnvs.txt')[sample_ID %in% samples$sample_ID, ]

examples <- rbind(dt[numsnp >= 40 & GT == 1 & prob >= 0.95, ][sample(1:.N, 2), ],
                  dt[numsnp >= 40 & GT == 2 & prob >= 0.95, ][sample(1:.N, 2), ],
                  dt[prob <= 0.05, ][sample(1:.N, 2), ])



# Run Grad CAM #

devtools::load_all()
# Load your trained model
model <- luz::luz_load('./joint.rds')

# Prepare PNGs
pngs_pt <- './tmp/gradcam_pngs'
snps <- fread('./data/hd_1kG_hg19.snppos.filtered.test.gz')
save_pngs_prediction(pngs_pt, examples, samples, snps, no_parall = F)

# Select a sample PNG
img_path <- list.files('./tmp/gradcam_pngs/', recursive = T, full.names = TRUE)[4]

# Generate and plot Grad-CAM
devtools::load_all()

gradcam_result <- generate_gradcam(model, img_path, target_class = 0)
plot_gradcam(gradcam_result)
gradcam_result <- generate_gradcam(model, img_path, target_class = 1)
plot_gradcam(gradcam_result)
gradcam_result <- generate_gradcam(model, img_path, target_class = 2)
plot_gradcam(gradcam_result)

# Or save to file
plot_gradcam(gradcam_result, save_path = './tmp/gradcam_visualization.png')
