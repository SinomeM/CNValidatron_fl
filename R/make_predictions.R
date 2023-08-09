
#' Save the PNG Plots of CNVs for Prediction
#'
#' This function uses a pre trained model to make a prediction on a new set
#' of CNVs
#'
#' @param model pre trained model loaded using `luz::luz_load()`
#' @param root root folder for the dataset, created using `save_pngs_prediction()`
#' @param cnvs cnv data.table in the usual format
#'
#' @export
#'
#' @import data.table
#' @import torchvision
#' @import luz
#' @import torch

save_pngs_prediction <- function(model, root, cnvs) {

  pred_dt <- image_folder_dataset(root, transform = . %>% transform_to_tensor())
  # the output is raw logits
  pred_tens <- predict(model, pred_dt)
  # convert to probabilities
  pred_tens <- 1 / (1 + exp(-pred_tens))

  # the classes are
  # false:   1
  # tru_del: 2
  # tru_dup: 3
  # unk_del: 4
  # unk_dup: 5


  pred_dt <- data.table(ix = uneval_dataset$samples[[1]],
                        class1 = round(as.numeric(uneval_pred[,1]), 5),
                        class2 = round(as.numeric(uneval_pred[,2]), 5),
                        class3 = round(as.numeric(uneval_pred[,3]), 5),
                        class4 = round(as.numeric(uneval_pred[,4]), 5),
                        class5 = round(as.numeric(uneval_pred[,5]), 5))

  for (i in 1:nrow(dt))
    dt[i, pred := which.max(c(class1, class2, class3, class4, class5))]

  dt[, sample_ID := gsub('.+new/', '', ix)]
  dt[, start := gsub('\\d+_', '', sample_ID)][, start := as.integer(gsub('.png', '', start))]
  dt[, sample_ID := as.integer(gsub('_\\w+.png', '', sample_ID))]
  dt[, ix := NULL]

  dt[, locus := paste0('locus', 1:nrow(dt))]
  dt <- merge(dt, cnvs[, .(sample_ID, chr, start, end, numsnp,
                           length, type, conf, batch, GT, CN)],
              by = c('sample_ID', 'start'))

}
