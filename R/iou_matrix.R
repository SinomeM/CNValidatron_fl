#' Compute the Intersection Over the Union for a set of CNVs
#'
#' Can be usefull for exploratory reasons but also as a base to
#' construct CNVRs.
#'
#' @param cnv usual CNVs `data.table`
#'
#' @import data.table
#'
#' @export

iou_matrix <- function(cnvs, chr_arms) {

  # index on center position
  dt_all <- copy(cnvs)
  dt_all[, center := start + length/2]

  # iterate over chromosomal arms to avoid having a too large vector
  out_list <- list()
  for (i in 1:nrow(chr_arms)) {

    a <- chr_arms[i]
    message(a$arm_ID, '\n')

    dt <- dt_all[chr == a$chr & start >= a$ start & end <= a$end, ]
    setorder(dt, center)
    dt[, cix := 1:.N]

    # create and fill the similarity matrix in the "long"
    # format, as a data.table. Better for modern R, e.g. ggplot2
    dt_s <- as.data.table(expand.grid(dt$cix, dt$cix))
    dt_s <- merge(dt_s, dt[, .(cix, start, end)], by.x = 'Var1', by.y = 'cix')
    dt_s <- merge(dt_s, dt[, .(cix, start, end)], by.x = 'Var2', by.y = 'cix')
    colnames(dt_s) <- c('cnvB', 'cnvA', 'stA', 'enA', 'stB', 'enB')

    dt_s[, iou := (pmin(enA, enB) - pmax(stA, stB)) /
                  (pmax(enA, enB) - pmin(stA, stB))]
    setorder(dt_s, cnvA)

    # output object is going to be a list
    out_list[[i]] <- dt_s[, .(cnvA, cnvB, iou, stA, stB, enA, enB)]
    gc()
  }
  names(out_list) <- chr_arms$arm_ID
  return(out_list)
}


#' IOU heatmap-like plot
#'
#' @param dt output of `iou_matrix()`, one chr arm at the time (one list item)
#'
#' @import data.table
#' @import ggplot2
#'
#' @export

plot_iou <- function(dt) {
   ggplot(dt, aes(cnvA, cnvB, fill = iou)) +
   scale_fill_gradient2(low = '#0072B2', high = '#D55E00', mid= 'white',
                        midpoint = 0, limit = c(-1,1), space = 'Lab') +
   geom_tile() + theme_minimal() + coord_fixed() +
   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
         axis.title.x=element_blank(), axis.text.y=element_blank(),
         axis.ticks.y=element_blank(), axis.title.y=element_blank())
}
