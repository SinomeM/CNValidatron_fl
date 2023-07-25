
# A function to lead the snps data from the tabix indexed intensity file.
# It will also perform some preporcessing, e.g. reduce the LRR interval to [-1.4 , 1.2],
# in hte future some denoising may be added

# @param cnv one line data.table in the usual cnv format
# @param samp one line samples file in the usual format
# @param snps the snps file for PennCNV, as data.table
# @param in_out_ratio ratio of bp otside the cnv vs inside, per side. A value of one means
#        the cnv length will be addded on each side

load_snps_tbx <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                          min_lrr, max_lrr, shrink_lrr = NULL) {
  chr <- cnv$chr
  start <- cnv$start
  end <- cnv$end
  len <- end - start + 1
  tbx_path <- samp$file_path_tabix

  st <- start - (in_out_ratio*len)
  st <- ifelse(st < 0, 0, st)

  dt <- fread(cmd = paste0("tabix ", tbx_path, " ", chr, ":", st,
                          "-", end + (in_out_ratio*len)), header = F)

  if (nrow(dt) == 0) stop('File: ', tbx_path, ' seems empty or broken\n')

  if (adjusted_lrr) colnames(dt) <- c("chr", "position", "end", "LRR", "BAF", "LRRadj")
  else colnames(dt) <- c("chr", "position", "LRR", "BAF")

  if (!is.null(snps)) dt <- dt[paste0(chr, position) %in% snps[, paste0(Chr, Position)], ]

  setorder(dt, position)

  if (adjusted_lrr) setnames(dt, c('BAF', 'LRRadj'), c('baf', 'lrr'))
  else setnames(dt, c('BAF', 'LRR'), c('baf', 'lrr'))

  ## some preprocessing ##
  # restrict the lrr space
  dt[lrr > max_lrr, lrr := max_lrr][lrr < min_lrr, lrr := min_lrr]
  # id lrr is missing outside of the cnv 'impute' it
  dt[is.na(lrr) & !between(position, start, end), lrr := dt[!between(position, start, end),
                                                              mean(lrr, na.rm = T)]]
  # if baf is missing outised the cnv set it to either 0 or 1
  dt[is.na(baf) & !between(position, start, end), baf := sample(rep(0:1, length.out = .N))]
  # if lrr or baf is missing inside the cnv exclude the point
  dt <- dt[!((is.na(lrr) | is.na(baf)) & between(position, start, end)),]

  ## some more processing ##
  dt[position < start, group := 1][between(position, start, end), group := 2][
       position > end, group := 3]
  ms1 <- dt[group == 1, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms2 <- dt[group == 2, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]
  ms3 <- dt[group == 3, c(mean(lrr, na.rm = T), sd(lrr, na.rm = T))]

  if (!is.null(shrink_lrr)) {
    # snps in each group get pulled towards the group mean relative to their distance
    dt[group == 1 & lrr > ms1[1], lrr := lrr - (lrr-ms1[1])*shrink_lrr][
         group == 1 & lrr < ms1[1], lrr := lrr + (lrr-ms1[1])*shrink_lrr]
    dt[group == 2 & lrr > ms2[1], lrr := lrr - (lrr-ms2[1])*shrink_lrr][
         group == 2 & lrr < ms2[1], lrr := lrr + (lrr-ms2[1])*shrink_lrr]
    dt[group == 3 & lrr > ms3[1], lrr := lrr - (lrr-ms3[1])*shrink_lrr][
         group == 3 & lrr < ms3[1], lrr := lrr + (lrr-ms3[1])*shrink_lrr]
  }

  # ouliers removal, 3SDs
  dt[group == 1 & !between(lrr, ms1[1]-3*ms1[2], ms1[1]+3*ms1[2]), lrr := NA]
  dt[group == 1 & !between(lrr, ms2[1]-3*ms2[2], ms2[1]+3*ms2[2]), lrr := NA]
  dt[group == 1 & !between(lrr, ms3[1]-3*ms3[2], ms3[1]+3*ms3[2]), lrr := NA]
  dt <- dt[!is.na(lrr), ]



  return(dt[, .(chr, position, lrr, baf)])

}

### --- ### --- ###

plot_cnv <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                     w = 64, z_ratio = 0.1, tmp_plot = 0, min_lrr = -1.2, max_lrr = 1,
                     shrink_lrr = NULL) {
  # initial checks
  if ((w %% 2) != 0) stop('w must be even')

  # compute w, z and k values
  z <- round(w * z_ratio)
  if ((z %% 2) != 0) z <- z - 1
  k <- (w - z)/2
  if (((k*2) + z) != w) stop('something wrong')

  # everything will be [0,w-1] then I will add 1 to make it [1,w]
  w <- w-1

  # load snps data
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)

  # move position to the x coordinates in the new system
  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len);  ee <- cnv$end + (in_out_ratio*len)
  dt[, x := round(((position-ss)/(ee-ss)) * w)]

  # each point need to be used for both lrr and baf so dt must be duplicated
  dt_lrr <- copy(dt)
  dt_baf <- copy(dt)

  # move lrr and baf on hte y coordinates in the new system
  dt_lrr[, y := round(((lrr-(min_lrr))/(max_lrr-(min_lrr))) * k)]
  dt_baf[, y := round(((baf-0)/(1-0)) * k) + k + z]

  # pixel coordinates must be > 0
  dt_baf[, ':=' (x = x+1, y = y+1)]
  dt_lrr[, ':=' (x = x+1, y = y+1)]
  w <- w+1

  # temporary check
  if (tmp_plot == 1) {
    a <- ggplot(dt_lrr, aes(x, y)) + geom_point() + xlim(0, w) + ylim(0, k) + theme_bw()
    b <- ggplot(dt_baf, aes(x, y)) + geom_point() + xlim(0, w) + ylim(k+z, w) + theme_bw()
    print(cowplot::plot_grid(b, a, ncol = 1))
  }

  # create the pixel map
  dt_lrr <- dt_lrr[, .N, by = c('x', 'y')]
  dt_baf <- dt_baf[, .N, by = c('x', 'y')]
  dt <- as.data.table(rbind(t(combn(1:w, 2)), t((combn(w:1, 2))))) # almost all combinations
  colnames(dt) <- c('x', 'y')
  dt <- rbind(dt, data.table(x = 1:w, y = 1:w)) # the diagonal was missing
  dt <- merge(dt, dt_lrr, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'a')
  dt <- merge(dt, dt_baf, by = c('x', 'y'), all.x = T)
  setnames(dt, 'N', 'b')
  dt[!is.na(a), value := a]
  dt[!is.na(b), value := b]
  dt[!is.na(a) & !is.na(b), value := a+b]
  dt[is.na(a) & is.na(b), value := 0]
  dt <- dt[, .(x, y, value)]

  if (tmp_plot == 2)
    print(ggplot(dt, aes(x, y, fill = value)) + geom_tile() + theme_bw() +
            scale_fill_gradient(low="white", high="black"))

  return(dt)

}

### --- ### --- ###

check_cnv <- function(cnv, samp, snps = NULL, in_out_ratio = 1, adjusted_lrr = T,
                      min_lrr = -1.2, max_lrr = 1, shrink_lrr = NULL) {

  # load snps data
  dt <- load_snps_tbx(cnv, samp, snps, in_out_ratio, adjusted_lrr, min_lrr, max_lrr, shrink_lrr)

  dt[between(position, cnv$start, cnv$end), inside := T][is.na(inside), inside := F]

  len <- cnv$end - cnv$start + 1
  ss <- cnv$start - (in_out_ratio*len)
  ee <- cnv$end + (in_out_ratio*len)

  a <- ggplot(dt, aes(x = position, y = baf, colour = inside)) + geom_point() +
    xlim(ss, ee) + theme_bw()
  b <- ggplot(dt, aes(x = position, y = lrr, colour = inside)) + geom_point() +
    xlim(ss, ee) + ylim(min_lrr, max_lrr) + theme_bw()

  cowplot::plot_grid(a,b, ncol = 1)

}


### --- ### --- ###

save_pngs_dataset <- function(root, cnvs, samps, snps, w = 64, in_out_ratio = 3,
                              shrink_lrr = 0.2, flip_chance = 0.5) {
  if (dir.exists(root)) stop('Root folder already exists. Delete existing folder or provide a dfferent path')

  dir.create(root)
  dir.create(paste0(root, '/true_del')); dir.create(paste0(root, '/true_dup'))
  dir.create(paste0(root, '/unk_dup')); dir.create(paste0(root, '/unk_del'))
  dir.create(paste0(root, '/false'))

  for (i in 1:nrow(cnvs)) {
    a <- cnvs[i]

    if (a$GT == 1)
      if (a$Visual_Output == 1) pt <- paste0(root, '/true_del/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 2) pt <- paste0(root, '/false/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 3) pt <- paste0(root, '/unk_del/', a$sample_ID, '_', a$start, '.png')
    if (a$GT == 2)
      if (a$Visual_Output == 1) pt <- paste0(root, '/true_dup/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 2) pt <- paste0(root, '/false/', a$sample_ID, '_', a$start, '.png')
      if (a$Visual_Output == 3) pt <- paste0(root, '/unk_dup/', a$sample_ID, '_', a$start, '.png')

      dt <- plot_cnv(a, samps[sample_ID == a[, sample_ID], ], snps = snps,
                     w = w, in_out_ratio = in_out_ratio, shrink_lrr = shrink_lrr)

      dt[, y := abs(y-(max(y)+1))] # to deal with how imager use the y axis
      imager::save.image(imager::as.cimg(dt), pt)

      # Data agumentation 1, image flipping
      if (runif(1) >= flip_chance) {
        dt[, x := abs(x-(max(x)+1))] # flip the x axis
        pt <- gsub('\\.png', '_flip\\.png', pt)
        imager::save.image(imager::as.cimg(dt), pt)
      }
  }
}
