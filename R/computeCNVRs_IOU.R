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
#' @param max_force_merge_rounds how many rounds of force merge are performed at the end of the function
#'
#' @import data.table
#' @import igraph
#'
#' @export

cnvrs_iou <- function(cnvs, chr_arm, screen_size = 500, min_iou = 0.75,
                      leiden_res = 1, plots_path = NA, min_n = 20,
                      max_force_merge_rounds = 5) {

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
        ig2 <- get_igraph_objs(dt2, min_iou, leiden_res, ii, cc$arm_ID, 'sml')
        dt2 <- ig2[[3]]
        dt_r2 <- create_cnvrs(dt2)
      }

      dt <- rbind(dt1, dt2); dt_r <- rbind(dt_r1, dt_r2)

      # Check for overlapping CNVRs and merge them
      out <- force_cnvr_merge(dt, dt_r)
      dt <- out[[1]]
      dt_r <- out[[2]]

      # update the output tables
      cnvs_with_CNVR <- rbind(cnvs_with_CNVR, dt)
      cnvrs <- rbind(cnvrs, dt_r)

      if (!is.na(plots_path)) save_igraph_plot(plots_path, ig[[1]], ig[[2]], ii)
      gc()
    }
  }

  for (i in 1:max_force_merge_rounds) {
    message('Final CNVRs merging round ', i, ' out of ', max_force_merge_rounds)
    dt <- force_cnvr_merge(cnvs_with_CNVR, cnvrs, verbose = T)
    cnvs_with_CNVR <- dt[[1]]
    cnvrs <- dt[[2]]
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
  setnames(dt_s, c('cix', 'start', 'end', 'i.cix', 'i.start', 'i.end'),
           c('cnvA', 'stA', 'enA', 'cnvB', 'stB', 'enB'))
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
  dt_g <- dt_g[, .(cnvA, cnvB, iou)]
  colnames(dt_g) <- c('a', 'b', 'weight')
  g <- graph_from_data_frame(dt_g, directed = F)
  gr <- cluster_leiden(g, resolution_parameter = leiden_res)
  a <- membership(gr)
  a <- data.table(cix = as.integer(names(a)), gr = a)

  # merge dt_r and igraph object here, careful, not all CNVs that got in are necessarily in dt_r
  dt_r <- merge(dt_r, a, by = 'cix')
  dt_r[, CNVR := paste0(arm, '_', ii, '_', ifelse(type != '', paste0(type, '_'), ''), gr)]
  dt_r[, gr := NULL]

  if (dt_m[, .N] > 0) {
    dt_m[, CNVR := paste0(arm, '_', ii, '_sngl_', 1:.N)]
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
    dt_r <- rbind(dt_r, data.table(CNVR = cr, chr = a[1, chr],
                                   start = a[, round(median(start))],
                                   end = a[, round(median(end))], n = a[, .N]))
  }
  return(dt_r)
}

force_cnvr_merge <- function(cnvs, cnvrs, verbose = F) {

  if (cnvs[CNVR %in% cnvrs[, CNVR], .N] != cnvs[, .N])
    stop('Some CNVs are assigned to missing CNVRs')
  if (cnvrs[CNVR %in% cnvs[, CNVR], .N] != cnvrs[, .N])
    stop('Some CNVRs are assigned to missing CNVs')

  cnvs[, length := end - start + 1]
  cnvrs[, length := end - start + 1]

  # run by chromosome
  for (i in cnvrs[, unique(chr)]) {
    if (verbose) message('Chr ', i)
    # cnvs and cnvrs per chr
    dt <- cnvs[chr == i, ]
    dt_r <- cnvrs[chr == i, ]
    setorder(dt_r, start)
    dt_r[, skip := F]
    # check each cnvr
    for (ii in 1:dt_r[, .N]) {
      a <- dt_r[ii]
      # skip if already merged
      if (a[, skip]) next

      # reciprocal overlap 0.75
      dt_r[, overlap := pmin(a$end, end) - pmax(a$start, start)]
      b <- dt_r[overlap >= 0.75*a$length & overlap >= 0.75*length, ]
      # skip if the only match is itself
      if (b[, .N] == 1) next

      dt_r[CNVR %in% b[, CNVR], skip := T]
      # update CNVs of to be merged CNVRs
      cnvs[CNVR %in% b[, CNVR], CNVR := a[, paste0(CNVR, '_fmrg')]]
    }
  }

  # recreate CNVRs from updated CNVs
  cnvrs <- CNValidatron:::create_cnvrs(cnvs)
  return(list(cnvs, cnvrs))
}

# not used at the moment
merge_cnvrs <- function(cnvrs, cnvs, min_iou, leiden_res, arm, ii) {

  # compute iou between all pairs of overlapping CNVRs
  setkey(cnvrs, start, end)
  dt <- foverlaps(cnvrs[, .(start, end, CNVR)], cnvrs[, .(start, end, CNVR)])
  setnames(dt, c('CNVR', 'start', 'end', 'i.CNVR', 'i.start', 'i.end'),
           c('cnvrA', 'stA', 'enA', 'cnvrB', 'stB', 'enB'))
  dt <- dt[cnvrA != cnvrB & !(is.na(stA) | is.na(stB)), ]
  dt[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
              (pmax(enA, enB) - pmin(stA, stB))]

  if (dt[, .N] == 0) return(list(cnvs, cnvrs))

  message('Merging come CNVRs!')

  # I can treat CNVRs as CNVs and use the same idea of the get_igraph_objs() function

  # cnvrs with connections
  dt_g <- dt[iou >= min_iou, ]
  setorder(dt_g, cnvrA)
  dt_r <- cnvrs[CNVR %in% dt_g[, unique(c(cnvrA, cnvrB))], ]
  dt_r[, CNVR_old := CNVR]
  # cnvs without connections
  dt_m <- cnvrs[!CNVR %in% dt_g[, unique(c(cnvrA, cnvrB))], ]

  # create the igraph network
  g <- graph_from_data_frame(dt_g[, .(cnvrA, cnvrB)], directed = F)
  gr <- cluster_leiden(g, resolution_parameter = leiden_res)
  a <- membership(gr)
  a <- data.table(CNVR = names(a), gr = a)
  # re construct the CNVRs table
  dt_r <- merge(dt_r, a, by = 'CNVR')
  dt_r[, CNVR := paste0(arm, '_', ii, '_mrg_', gr)]
  dt_r[, gr := NULL]

  # Updated CNVRs info in the CNVs table
  for (i in dt_r[, unique(CNVR)])
    cnvs[CNVR %in% dt_r[CNVR == i, CNVR_old], CNVR := i]

  # recreate CNVRs table is from updated CNVs table
  cnvrs <- create_cnvrs(cnvs)

  return(list(cnvs, cnvrs))
}
