
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

  dt_rbind <- data.table()

  for (bb in 1:batches) {
    message(paste0('Batch ', bb, ' of ', batches))
    pt <- paste0(root, '/batch', bb, '/')

    pred_dt <- image_folder_dataset(pt, transform = . %>% transform_to_tensor())
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
    colnames(pred_probs) <- c('p_false', 'p_true_del', 'p_true_dup', 'p_unk_del', 'p_unk_dup')
    pred_dt <- cbind(pred_dt, pred_probs)

    dt_rbind <- rbind(dt_rbind, pred_dt)
  }


  dt_rbind[, sample_ID := gsub('samp\\w+', '', ix)][,
            sample_ID := as.character(gsub('_.+', '', sample_ID))][,
            start := gsub('_st\\d+', '', ix)][,
            start := as.integer(gsub('_st', '', start))][,
            real_numsnp := gsub('nsnp\\w+', '', ix)][,
            real_numsnp := as.integer(gsub('nsnp', '', real_numsnp))]

  if (return_pred_dt)
    returt(dt_rbind)

  dt_rbind[, ix := NULL]
  cnvs[, ':=' (sample_ID = as.character(sample_ID), start = as.integer(start))]

  dt_rbind <- merge(dt_rbind, cnvs[, .(sample_ID, chr, start, end, numsnp,
                                length, conf, GT, CN)],
                   by = c('sample_ID', 'start'))

  returt(dt_rbind)
}
