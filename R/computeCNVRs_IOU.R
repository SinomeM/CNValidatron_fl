#' Create CNVRs in a set of CNVS based on network analysis
#'
#' A network is defined a set of CNVs that overlap with each other. Usually there
#' are several ones in each chromosomal arm. Given a minimum IOU value
#' each pairs of CNVs is considered as two connected nodes. The CNVRs are
#' computed using community the detection algorithm 'leiden'
#'
#' @param cnvs usual CNVs `data.table`
#' @param chr_arms chromsome arms location, from `QCtreeCNV` package
#' @param screen_size size of the comb that screen the cnvs to look for individual networks
#' @param min_iou minimum IOU filter. Define how similar two CNVs must be in order to be considered connected
#' @param leiden_res lower value will tend to create more smaller communities
#' @param plot_path not yet implemented
#' @param min_n min n for recomputing CNVRs (small CNVRs might "drown" in very large networks)
#'
#' @import data.table
#' @import igraph
#'
#' @export

cnvrs_iou <- function(cnvs, chr_arm, screen_size = 500, min_iou = 0.75,
                      leiden_res = 1, plots_path = NA, min_n = 20) {

  cnvs[, center := round(start + (end-start+1)/2)]
  cnvs_with_CNVR <- data.table()
  cnvrs <- data.table()

  for (i in 1:nrow(chr_arm)) {

    cc <- chr_arm[i]
    message(cc$arm_ID)
    cnvs_arm <- cnvs[chr == cc$chr & (between(start, cc$start, cc$end) | between(end, cc$start, cc$end)), ]

    if (cnvs_arm[, .N] == 0) next

    splits <- create_splits_foverlaps(cnvs_arm, cc, screen_size)

    # now we have a more manageable set of CNVs we can proceed with the IOU matrix and the network analysis
    for (ii in 1:nrow(splits)) {
      a <- splits[ii]
      dt <- cnvs_arm[start >= a$st & end <= a$en, ]
      if (nrow(dt) == 0) next
      message('Analysing network ', ii, ', ', nrow(dt), ' CNVs...')

      # create igraph object and groups
      ig <- get_igraph_objs(dt, min_iou, leiden_res, ii, cc$arm_ID)

      # reconstruct CNVRs
      dt <- ig[[3]]
      dt_r <- create_cnvrs(dt)

      # reprocess all CNVs from CNVRs with n < 10
      dt_r1 <- dt_r[n >= min_n, ]
      dt_r2 <- dt_r[n < min_n, ]
      dt1 <- dt[CNVR %in% dt_r1$CNVR, ]
      dt2 <- dt[CNVR %in% dt_r2$CNVR, ]

      if (nrow(dt2) > 1) {
        ig2 <- get_igraph_objs(dt2, min_iou, leiden_res, ii, cc$arm_ID, 'small')
        # dt2[, CNVR := paste0(cc$arm_ID, '_', ii, 'small_', membership(ig2[[2]]))]
        dt2 <- ig2[[3]]
        dt_r2 <- create_cnvrs(dt2)
      }

      dt <- rbind(dt1, dt2); dt_r <- rbind(dt_r1, dt_r2)

      # Check for overlapping CNVRs and merge them
      while(T) {
        # TO BE TESTED !!!!
        n <- nrow(dt_r)
        out <- merge_cnvrs(dt_r, dt, min_iou)
        dt <- out[[1]]
        dt_r <- out[[2]]
        n1 <- nrow(dt_r)
        if(n1 < n) message('merged some CNVRs')
        else break
      }

      # update the output tables
      cnvs_with_CNVR <- rbind(cnvs_with_CNVR, dt)
      cnvrs <- rbind(cnvrs, dt_r)

      if (!is.na(plots_path)) save_igraph_plot(plots_path, ig[[1]], ig[[2]], ii)
      gc()
    }
  }

  return(list(cnvs_with_CNVR, cnvrs))
}


create_splits <- function(cnvs_arm, window_size, cc) {
  # create the bins
  breaks <- seq(from = cnvs_arm[, min(start)], to = cnvs_arm[, max(end)], by = window_size)
  breaks[length(breaks)] <- cnvs_arm[, max(end) + 1]
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

create_splits_foverlaps <- function(cnvs_arm, cc, screen_size = 500) {

  # create a sort of "comb" with the teeth at screen_size distance
  breaks <- seq(from = cc$start, to = cc$end, by = screen_size)
  bins <- data.table(start = breaks, end = breaks + 1, ix = 1:length(breaks))

  # foverlaps between the comb and the CNVs
  setkey(cnvs_arm, start, end)
  dt <- foverlaps(bins, cnvs_arm[, .(start, end)])
  # the teeth with no match (NA in start and end) can be a split boundary
  dt <- dt[is.na(start), ]
  # keep only those were at least one is skipped, meaning the previous did overlap wit a CNV
  dt[, ixp := c(0, dt$ix[-length(dt$ix)])]
  splits <- unique(rbind(dt[1, .(i.start, ix)], dt[ixp < ix - 1, .(i.start, ix)], dt[.N, .(i.start, ix)]))
  colnames(splits) <- c('start', 'ix')
  # start of of following is the end of previous, this result in .N-1 splits
  splits <- splits[1:(.N-1),][, end := (splits$start[-1]) - 1]
  message(nrow(splits), ' different possible networks detected in ', nrow(cnvs_arm), ' CNVs')

  return(splits)
}

get_igraph_objs <- function(dt, min_iou, leiden_res, ii, arm, type = '') {

  # get all overlapping CNVs
  setkey(dt, start, end)
  dt[, cix := 1:.N]
  dt_s <- foverlaps(dt[, .(cix, start, end)], dt[, .(cix, start, end)])
  colnames(dt_s) <- c('cnvA', 'stA', 'enA', 'cnvB', 'stB', 'enB')
  dt_s <- dt_s[cnvA != cnvB, ]

  dt_s[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
                (pmax(enA, enB) - pmin(stA, stB))]

  # cnvs with connections
  dt_g <- dt_s[iou >= min_iou, ]
  setorder(dt_g, cnvA)
  dt_r <- dt[cix %in% dt_g[, unique(c(cnvA, cnvB))], ]
  # cnvs without connections
  dt_m <- dt[!cix %in% dt_g[, unique(c(cnvA, cnvB))], ]

  # create the igraph network
  g <- graph_from_data_frame(dt_g[, .(cnvA, cnvB)], directed = F)
  gr <- cluster_leiden(g, resolution_parameter = leiden_res)
  a <- membership(gr)
  a <- data.table(cix = as.integer(names(a)), gr = a)

  # merge dt_r and igraph object here, careful, not all CNVs that got in are necessarily in dt_r
  dt_r <- merge(dt_r, a, by = 'cix')
  dt_r[, CNVR := paste0(arm, '_', ii, '_', ifelse(type != '', paste0(type, '_'), ''), gr)]
  dt_r[, gr := NULL]

  if (dt_m[, .N] > 0) {
    dt_m[, CNVR := paste0(arm, '_', ii, '_singleton_', 1:.N)]
    dt <- rbind(dt_m, dt_r)
    return(list(g, gr, dt))
  }

  return(list(g, gr, dt_r))
}

save_igraph_plot <- function(plots_path, g, gr, ii) {
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
  continue <- T; i <- 0; max <- round(nrow(cnvrs)/5)
  cnvrs[, len := end - start]
  while(i <= max) {
    # select one random cnvr, from the largest half
    a <- cnvrs[len >= cnvrs[, mean(len)/2], ][sample(1)]
    len <- a[, end - start + 1]
    # check if there is any other CNVR with iou >= min_iou
    cnvrs[, iou := (pmin(a$end, end) - pmax(a$start, start)) /
                   (pmax(a$end, end) - pmin(a$start, start))]
    b <- cnvrs[iou >= min_iou, ]

    if (nrow(b) > 1) {
      # merge all CNVRs and update relative CNVs
      c <- fsetdiff(cnvrs, b)
      b[1][, ':=' (start = min(start), end = max(end), n = sum(n))]
      cnvs[CNVR %chin% b$CNVR, CNVR := b[1, paste0(CNVR, '_merge')], ]
      cnvrs <- rbind(c, b[1])
    }

    i <- i + 1
  }
  return(list(cnvs, cnvrs))
}
