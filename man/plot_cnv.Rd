% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_cnv_pixel_image.R
\name{plot_cnv}
\alias{plot_cnv}
\title{Create the matrix for the PNG image}
\usage{
plot_cnv(
  cnv,
  samp,
  snps = NULL,
  adjusted_lrr = T,
  tmp_plot = 0,
  min_lrr = -1.4,
  max_lrr = 1.3,
  simple_min_max = F,
  use_log = F,
  blank_small = F,
  shrink_lrr = 0.1,
  w = 96,
  z = 4,
  k1 = 31,
  k2 = 26,
  l_wind = 2e+07,
  mx_lr = 2
)
}
\arguments{
\item{cnv}{see load_snps_tbx() documentation}

\item{samp}{see load_snps_tbx() documentation}

\item{snps}{see load_snps_tbx() documentation}

\item{adjusted_lrr}{see load_snps_tbx() documentation}

\item{tmp_plot}{for developing, if set to 1 plot the "normal"
LRR/BAF plot, if set to 2 plot the pixelated image in R}

\item{min_lrr}{see load_snps_tbx() documentation}

\item{max_lrr}{see load_snps_tbx() documentation}

\item{shrink_lrr}{see load_snps_tbx() documentation}
}
\description{
This function create the pixel matrix that can be saved
as a PNG for further use. All processing is done in the
function load_snps_tbx().
}
