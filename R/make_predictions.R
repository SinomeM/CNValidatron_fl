
#' Make Prediction on a set of CNVs
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

make_predicitons <- function(model, root, cnvs) {

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


  pred_dt <- data.table(ix = pred_dt$samples[[1]],
                        class1 = round(as.numeric(pred_tens[,1]), 5),
                        class2 = round(as.numeric(pred_tens[,2]), 5),
                        class3 = round(as.numeric(pred_tens[,3]), 5),
                        class4 = round(as.numeric(pred_tens[,4]), 5),
                        class5 = round(as.numeric(pred_tens[,5]), 5))

  for (i in 1:nrow(pred_dt))
    pred_dt[i, pred := which.max(c(class1, class2, class3, class4, class5))]

  pred_dt[, sample_ID := gsub('.+new/', '', ix)][,
            start := gsub('\\d+_', '', sample_ID)][,
            start := as.integer(gsub('.png', '', start))][,
            sample_ID := as.integer(gsub('_\\w+.png', '', sample_ID))][,
            ix := NULL]

  pred_dt <- merge(pred_dt, cnvs[, .(sample_ID, chr, start, end, numsnp,
                                length, type, conf, batch, GT, CN)],
                   by = c('sample_ID', 'start'))

}
