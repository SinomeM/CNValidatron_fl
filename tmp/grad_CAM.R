

# Grad CAM plots on some examples #

# Find very good CNVs to use as examples, both deletions and duplications
files <- gsub('\\.tabix.+', '', list.files('./data/1000Genomes_data/tbx/'))
samples <- fread('./data/1000Genomes_data/common_samples.txt')[hd_id %in% files, ]
samples <- samples[, .(sample_ID = Sample_Group,
            file_path_tabix = paste0('/Users/sm/Documents/GitHub/CNValidatron_fl/data/1000Genomes_data/tbx/',
            hd_id, '.tabix.gz'))]
dt <- fread('./data/1000Genomes_data/1kG_cnvs.txt')[sample_ID %in% samples$sample_ID, ]

good_del <- dt[GT == 1 & prob >= 0.95, ][sample(1:.N, 1), ]
good_dup <- dt[GT == 2 & prob >= 0.95, ][sample(1:.N, 1), ]
bad_cnv <- dt[prob <= 0.05, ][sample(1:.N, 1), ]


# 
devtools::load_all()