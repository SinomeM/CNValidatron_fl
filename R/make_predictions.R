
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

make_predictions <- function(model, root, cnvs, return_pred_dt = F, batches = 1000) {


  pred_dt <- image_folder_dataset(root, transform = . %>% transform_to_tensor())
  # the output is raw logits
  pred_tens <- predict(model, pred_dt)
  # convert to probabilities
  pred_tens <- nnf_softmax(pred_tens, dim = 2)

  # the classes are
  # false:   1
  # tru_del: 3
  # tru_dup: 3

  pred_ix <- as.matrix(torch_argmax(pred_tens, dim = 2))[, 1]
  pred_probs <- round(as.matrix(pred_tens), 3)

  pred_dt <- data.table(ix = pred_dt$samples[[1]], pred = pred_ix)

  for (i in 1:nrow(pred_dt))
    pred_dt[i, pred_prob := pred_probs[i, pred_ix[i]]]

  pred_probs <- as.data.table(pred_probs)
  colnames(pred_probs) <- c('p_false', 'p_true_del', 'p_true_dup')
  pred_dt <- cbind(pred_dt, pred_probs)


  pred_dt[, sample_ID := gsub('.+samp', '', ix)][,
            sample_ID := as.character(gsub('_.+', '', sample_ID))][,
            start := gsub('.+st', '', ix)][,
            start := as.integer(gsub('_.+', '', start))][,
            real_numsnp := gsub('.+nsnp', '', ix)][,
            real_numsnp := as.integer(gsub('\\.png', '', real_numsnp))]

  if (return_pred_dt)
    returt(pred_dt)

  pred_dt[, ix := NULL]
  cnvs[, ':=' (sample_ID = as.character(sample_ID), start = as.integer(start))]

  pred_dt <- merge(pred_dt, cnvs[, .(sample_ID, chr, start, end, numsnp,
                                length, conf, GT, CN)],
                   by = c('sample_ID', 'start'))

  returt(pred_dt)
}
