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
#' @import igraph
#'
#' @export

cnvrs_iou <- function(cnvs, chr_arm, window_size = 500000, min_iou = 0.75,
                      leiden_res = 1, plots_path = NA, min_n = 10) {

  cnvs[, center := start + (end-start+1)/2]
  cnvs_with_CNVR <- data.table()
  cnvrs <- data.table()

  for (i in 1:nrow(chr_arm)) {

    cc <- chr_arm[i]
    message(cc$arm_ID)
    cnvs_arm <- cnvs[chr == cc$chr & (between(start, cc$start, cc$end) | between(end, cc$start, cc$end)), ]

    if (cnvs_arm[, .N] == 0) next

    splits <- create_splits(cnvs_arm, window_size, cc)

    # now we have a more manageable set of CNVs we can proceed with the IOU matrix and the network analysis
    for (ii in 1:nrow(splits)) {
      a <- splits[ii]
      dt <- cnvs_arm[start >= a$st & end <= a$en, ]
      if (nrow(dt) == 0) next
      message('Analysing network ', ii, ', ', nrow(dt), ' CNVs...')

      # create igraph object and groups
      ig <- get_igraph_objs(dt, min_iou, leiden_res)
      # assign each CNV to a CNVR
      dt[, CNVR := paste0(cc$arm_ID, '_', ii, '_', membership(ig[[2]]))]

      # reconstruct CNVRs
      dt_r <- create_cnvrs(dt)

      # reprocess all CNVs from CNVRs with n < 10
      dt_r1 <- dt_r[n >= min_n, ]
      dt_r2 <- dt_r[n < min_n, ]
      dt1 <- dt[CNVR %in% dt_r1$CNVR, ]
      dt2 <- dt[CNVR %in% dt_r2$CNVR, ]

      if (nrow(dt2) > 1) {
        ig2 <- get_igraph_objs(dt2, min_iou, leiden_res)
        dt2[, CNVR := paste0(cc$arm_ID, '_', ii, 'small_', membership(ig2[[2]]))]
        dt_r2 <- create_cnvrs(dt2)
      }

      dt <- rbind(dt1, dt2); dt_r <- rbind(dt_r1, dt_r2)
      # Check for overlapping CNVRs and merge them
      n <- 2; n1 <- 1 # initialise the loop
      while(n1 < n) {
        # TO BE TESTED !!!!
        n <- nrow(dt_r)
        out <- merge_cnvrs(dt, dt_r, min_iou)
        dt <- out[[1]]
        dt_r <- out[[2]]
        n1 <- nrow(dt_r)
      }

      # update the output tables
      cnvs_with_CNVR <- rbind(cnvs_with_CNVR, dt)
      cnvrs <- rbind(cnvrs, dt_r)

      if (!is.na(plots_path)) save_igraph_plot(ig[[1]], ig[[2]])
      gc()
    }
  }


  # get back the CNVs that are excluded in the CNVRs detection (no connections iou > min_iou)
  # ... #

  return(list(cnvs_with_CNVR, cnvrs))
}


create_splits <- function(cnvs_arm, window_size, cc) {
  # create the bins
  breaks <- seq(from = cnvs_arm[, min(start)], to = cnvs_arm[, max(end)], by = window_size)
  bins <- data.table(st = breaks[-length(breaks)], en = breaks[-1], ix = 1:(length(breaks)-1))

  # count CNVs per bins
  for (ii in 1:nrow(bins))
    bins[ii, N := cnvs_arm[chr == cc$chr & (between(start, st, en) | between(end, st, en)), .N] ]

  # get 0 count bins and create the splits
  sp <- c(cc$start, bins[N == 0, st + (en-st+1)/2], cc$end)
  splits <- data.table(st = sp[-length(sp)], en = sp[-1], ix = 1:(length(sp)-1))
  message(nrow(splits), ' different possible networks detected in ', nrow(cnvs_arm), ' CNVs')

  return(splits)
}

get_igraph_objs <- function(dt, min_iou, leiden_res) {
  setorder(dt, center)
  dt[, cix := 1:.N]

  # create and fill the similarity matrix in the "long" format
  dt_s <- as.data.table(expand.grid(dt$cix, dt$cix))
  dt_s <- merge(dt_s, dt[, .(cix, start, end)], by.x = 'Var1', by.y = 'cix')
  dt_s <- merge(dt_s, dt[, .(cix, start, end)], by.x = 'Var2', by.y = 'cix')
  colnames(dt_s) <- c('cnvB', 'cnvA', 'stA', 'enA', 'stB', 'enB')

  dt_s[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
                (pmax(enA, enB) - pmin(stA, stB))]

  dt_s <- dt_s[iou >= min_iou, ]
  setorder(dt_s, cnvA)

  # create the igraph network
  g <- graph_from_data_frame(dt_s[,1:2], directed = F)
  gr <- cluster_leiden(g, resolution_parameter = leiden_res)

  return(list(g, gr))
}

save_igraph_plot <- function(g, gr) {
  colors <- rainbow(max(membership(gr)))
  pl <- plot(g, vertex.color = colors[membership(gr)],
             layout = layout_nicely, vertex.size = 5, vertex.label.cex = 0.3, arrow.mode = 0)
  ggsave(pl, paste0(plots_path, '/', cc$arm_ID, '_', ii, '.png'))
}

create_cnvrs <- function(dt) {
  dt_r <- data.table()
  for (cr in dt[, unique(CNVR)]) {
    a <- dt[CNVR == cr, ]
    dt_r <- rbind(dt_r, data.table(CNVR = cr, chr = a[1, chr], start = a[, min(start)], end = a[, max(end)], n = a[, .N]))
  }
  return(dt_r)
}

# TO BE TESTED!!!!
merge_cnvrs <- function(cnvrs, cnvs, min_iou) {
  # initialise while
  continue == T; i <- 1; max <- nrow(cnvr)
  while(continue) {
    # select one random cnvr
    a <- cnvrs[sample(1)]
    len <- a[, end - start + 1]
    # check if there is any other CNVR with iou >= min_iou
    cnvrs[, iou := (pmin(a$end, end) - pmax(a$start, start)) /
                   (pmax(a$end, end) - pmin(a$start, start))]
    b <- cnvrs[iou >= min_iou, ]

    if (nrow(b) > 1) {
      # merge all CNVRs and update relative CNVs
      c <- fsetdiff(cnvrs, b)
      b[1][, ':=' (start = min(start), end = max(end), n = sum(n))]
      cnvs[CNVR %chin% b$CNVR, CNVR := b[1, CNVR], ]
    }

    # exit the while after a certain number of loops
    if (i > max) break
    i <- 1 + 1
  }
  return(list(cnvs, cnvrs))
}
