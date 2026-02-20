
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

make_predictions <- function(model, root, cnvs, return_pred_dt = F,
                             clean_out = T) {

  pred_dt_rbind <- data.table()
  samples <- unique(cnvs$sample_ID)

  for (i in samples) {

    bpt <- paste0(root, '/', i, '/')
    # if the folder does not exist, or it is empty, skip
    if (!dir.exists(bpt)) next
    if (length(list.files(bpt)) == 0) next

    pred_dt <- image_folder_dataset(bpt, transform = . %>% transform_to_tensor())
    # the output is raw logits
    pred_tens <- predict(model, pred_dt)

    # normal case, multiple CNVs per samples
    if (length(dim(pred_tens)) == 2) {
      # convert to probabilities
      pred_tens <- nnf_softmax(pred_tens, dim = 2)

      # the classes are
      # false:   1
      # tru_del: 2
      # tru_dup: 3
      pred_ix <- as.matrix(torch_argmax(pred_tens, dim = 2))[, 1]
      pred_probs <- round(as.matrix(pred_tens), 3)

      pred_dt <- data.table(ix = pred_dt$samples[[1]], pred = pred_ix)

      for (i in 1:nrow(pred_dt))
        pred_dt[i, pred_prob := pred_probs[i, pred_ix[i]]]

      pred_probs <- as.data.table(pred_probs)
      colnames(pred_probs) <- c('p_false', 'p_true_del', 'p_true_dup')
    }

    # edge case, only one CNV for the sample
    if (length(dim(pred_tens)) == 1) {
      pred_tens <- nnf_softmax(pred_tens, dim = 1)
      pred_ix <- as.matrix(torch_argmax(pred_tens, dim = 1))[, 1]
      pred_probs <- round(as.matrix(pred_tens), 3)
      pred_dt <- data.table(ix = pred_dt$samples[[1]], pred = pred_ix)

      for (i in 1:nrow(pred_dt))
        pred_dt[i, pred_prob := pred_probs[pred_ix[i]]]

      pred_probs <- data.table(p_false = pred_probs[1],
                               p_true_del = pred_probs[2],
                               p_true_dup = pred_probs[3])
    }
    
    pred_dt <- cbind(pred_dt, pred_probs)

    pred_dt_rbind <- rbind(pred_dt_rbind, pred_dt)
  }

  if (return_pred_dt)
    return(pred_dt_rbind)

  pred_dt_rbind[, ':=' (sample_ID = as.character(gsub(".*/([^/]+)/new/.*", "\\1", ix)),
                        chr = as.integer(gsub(".*chr([0-9]+)_.*", "\\1", ix)),
                        start = as.integer(gsub(".*_st([0-9]+)_.*", "\\1", ix)),
                        real_numsnp = as.integer(gsub(".*_nsnp([0-9]+)\\.png", "\\1", ix)))]

  pred_dt_rbind[, ix := NULL]
  cnvs[, ':=' (sample_ID = as.character(sample_ID),
               chr = as.integer(chr),
               start = as.integer(start))]

  pred_dt <- merge(pred_dt_rbind, cnvs,
                   by = c('sample_ID', 'chr', 'start'))

  # Assign the probabilty based on the CNV type (GT)
  pred_dt[GT == 1, prob := p_true_del]
  pred_dt[GT == 2, prob := p_true_dup]

  if (clean_out)
    pred_dt[, c('p_false', 'p_true_del', 'p_true_dup', 'pred', 'pred_prob') := NULL]

  return(pred_dt)
}
