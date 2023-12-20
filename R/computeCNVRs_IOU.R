#' Create CNVRs in a set of CNVS based on network analysis
#'
#' A network is defined a set of CNVs that overlap with each other. Usually there
#' are several ones in each chromosomal arm. Given a minimum IOU value
#' each pairs of CNVs is considered as two connected nodes. The CNVRs are
#' computed using community the detection algorithm 'leiden'
#'
#' @param cnv usual CNVs `data.table`
#' @param chr_arms chromsome arms location, from `QCtreeCNV` package
#' @param min_iou minimum IOU filter. Define how similar two CNVs must be in order to be considered connected
#' @param leiden_res lower value will tend to create more smaller communities
#' @param plot_path not yet implemented
#'
#' @import data.table
#'
#' @export

cnvrs_iou <- function(cnvs, chr_arm, window_size = 500000, min_iou = 0.75, leiden_res = 1, plots_path = NA) {

  cnvs[, center := start + (end-start+1)/2]
  cnvs_with_CNVR <- data.table()

  for (i in 1:nrow(chr_arm)) {

    cc <- chr_arm[i]
    message(cc$arm_ID)
    cnvs_arm <- cnvs[chr == cc$chr & (between(start, cc$start, cc$end) | between(end, cc$start, cc$end)), ]

    if (cnvs_arm[, .N] == 0) next

    # create the bins
    breaks <- seq(from = cnvs_arm[, min(start)], to = cnvs_arm[, max(end)], by = window_size)
    bins <- data.table(st = breaks[-length(breaks)], en = breaks[-1], ix = 1:(length(breaks)-1))

    # count CNVs per bins
    for (ii in 1:nrow(bins))
      bins[ii, N := cnvs[chr == cc$chr & (between(start, st, en) | between(end, st, en)), .N] ]

    # get 0 count bins and create the splits
    sp <- c(cc$start, bins[N == 0, st + (en-st+1)/2], cc$end)
    splits <- data.table(st = sp[-length(sp)], en = sp[-1], ix = 1:(length(sp)-1))
    message(nrow(splits), ' different networks detected in ', nrow(cnvs_arm), ' CNVs')

    # now we have a more manageable set of CNVs we can proceed with the IOU matrix and the network analysis
    for (ii in 1:nrow(splits)) {
      a <- splits[ii]
      dt <- cnvs_arm[start >= a$st & end <= a$en, ]
      message('Analysing network ', ii, ', ', nrow(dt), ' CNVs...')

      setorder(dt, center); dt[, cix := 1:.N]

      # create and fill the similarity matrix in the "long" format
      dt_s <- as.data.table(expand.grid(dt$cix, dt$cix))
      dt_s <- merge(dt_s, dt[, .(cix, start, end)], by.x = 'Var1', by.y = 'cix')
      dt_s <- merge(dt_s, dt[, .(cix, start, end)], by.x = 'Var2', by.y = 'cix')
      colnames(dt_s) <- c('cnvB', 'cnvA', 'stA', 'enA', 'stB', 'enB')

      dt_s[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
                    (pmax(enA, enB) - pmin(stA, stB))]

      dt_s <- dt_s[iou >= min_iou, ]; setorder(dt_s, cnvA)

      # create the igraph network
      g <- graph_from_data_frame(dt_s[,1:2], directed = F)
      gr <- cluster_leiden(g, resolution_parameter = leiden_res)

      # add option to save the plot
      if (!is.na(plots_path)) {
        colors <- rainbow(max(membership(gr)))
        pl <- plot(g, vertex.color = colors[membership(gr)],
                   layout = layout_nicely, vertex.size = 5, vertex.label.cex = 0.3, arrow.mode = 0)
        ggsave(pl, paste0(cc$arm_ID, '_', ii, '.png'))
      }

      dt[, CNVR := paste0(cc$arm_ID, '_', ii, '_', membership(gr))]
      cnvs_with_CNVR <- rbind(cnvs_with_CNVR, dt)

      gc()
    }
  }
  return(cnvs_with_CNVR)
}
