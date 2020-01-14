
#' Run `eBayes` on NIfTI Images
#'
#' @inheritParams nifti_lmFit
#' @param proportion numeric value between 0 and 1, assumed proportion of
#' voxels which are differentially expressed, passed to \code{\link{eBayes}}
#' @param stdev.coef.lim numeric vector of length 2, assumed
#' lower and upper limits for the standard deviation of
#' log2-fold-changes for differentially expressed voxels,
#' passed to \code{\link{eBayes}}
#' @param trend ogical, should an intensity-trend be allowed for the
#' prior variance? Default is that the prior variance is constant,
#' passed to \code{\link{eBayes}}
#' @param robust logical, should the estimation of df.prior and
#' var.prior be robustified against outlier sample variances,
#' passed to \code{\link{eBayes}}
#' @param winsor.tail.p numeric vector of length 1 or 2,
#' giving left and right tail proportions of x to Winsorize.
#' Used only when \code{robust=TRUE}, passed to \code{\link{eBayes}}
#' @param adjust.method method to adjust p-value,
#'  passed to \code{\link{topTable}}
#' @param coef column number or column name specifying which
#' coefficient or contrast of the linear model is of interest,
#' passed to \code{\link{topTable}}
#' @return
#' @export
#'
#' @return A list of output images from `eBayes` and the linear fit
#' @examples
nifti_eBayes = function(
  imgs, mask, verbose = TRUE,
  ...,
  proportion = 0.01, stdev.coef.lim = c(0.1,4),
  trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1),
  coef = NULL,
  adjust.method = "BH") {

  img = RNifti::asNifti(imgs[[1]])
  if (!is.null(mask)) {
    mask =  RNifti::asNifti(mask) > 0
  } else {
    mask = array(TRUE, dim = dim(img))
  }

  mat = nifti_images_to_matrix(imgs, mask, verbose = verbose)
  n = nrow(mat)
  fit = limma::lmFit(mat, ...)

  eb.fit = limma::eBayes(
    fit,
    proportion = proportion,
    stdev.coef.lim = stdev.coef.lim,
    trend = trend, robust = robust,
    winsor.tail.p = winsor.tail.p)

  cols_to_grab = c("coefficients", "stdev.unscaled",
                   "t", "p.value", "lods", "adjusted_p_value")
  eb.fit$adjusted_p_value = limma::topTable(
    eb.fit,
    sort.by = "none",
    coef = coef,
    adjust.method = adjust.method,
    number = n)$adj.P.Val
  out_images = lapply(cols_to_grab, function(x) {
    remake_nifti_image(eb.fit[[x]], img = img, mask = mask)
  })
  names(out_images) = cols_to_grab

  L = list(
    images = out_images,
    empirical_bayes = eb.fit,
    adjust.method = adjust.method,
    lm_fit = fit)
  L
}
