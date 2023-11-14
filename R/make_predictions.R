
#' Make Prediction on a set of CNVs
#'
#' This function uses a pre trained model to make a prediction on a new set
#' of CNVs
#'
#' @param model pre trained model loaded using `luz::luz_load()`
#' @param root root folder for the dataset, created using `save_pngs_prediction()`
#' @param cnvs cnv data.table in the usual format
#' @param return_pred_dt returt the prediction table before merging w/ the cnv table, set
#'          it to T if there is a problem with the `sample_ID`
#'
#' @export
#'
#' @import data.table
#' @import torchvision
#' @import luz
#' @import torch

make_predictions <- function(model, root, cnvs, return_pred_dt = F) {

  pred_dt <- image_folder_dataset(root, transform = . %>% transform_to_tensor())
  # the output is raw logits
  pred_tens <- predict(model, pred_dt)
  # convert to probabilities
  pred_tens <- nnf_softmax(pred_tens, dim = 2)

  # the classes are
  # false:   1
  # tru_del: 2
  # tru_dup: 3
  # unk_del: 4
  # unk_dup: 5

  pred_ix <- as.matrix(torch_argmax(pred_tens, dim = 2))[, 1]
  pred_probs <- round(as.matrix(pred_tens), 3)

  pred_dt <- data.table(ix = pred_dt$samples[[1]], pred = pred_ix)

  for (i in 1:nrow(pred_dt))
    pred_dt[i, pred_prob := pred_probs[i, pred_ix[i]]]

  pred_probs <- as.data.table(pred_probs)
  colnames(pred_probs) <- paste0('prob', 1:5)
  pred_dt <- cbind(pred_dt, pred_probs)

  pred_dt[, sample_ID := gsub('.+new/', '', ix)][,
            start := gsub('\\d+_', '', sample_ID)][,
            start := as.integer(gsub('.png', '', start))][,
            sample_ID := as.integer(gsub('_\\w+.png', '', sample_ID))][,
            ix := NULL]

  pred_dt <- merge(pred_dt, cnvs[, .(sample_ID, chr, start, end, numsnp,
                                length, conf, GT, CN)],
                   by = c('sample_ID', 'start'))

}
