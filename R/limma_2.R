library(RNifti)
library(limma)


limma_images = function(
  imgs, mask, verbose = TRUE,
  ...,
  proportion = 0.01, stdev.coef.lim = c(0.1,4),
  trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1),
  adjust.method = "BH") {

  img = asNifti(imgs[[1]])
  if (!is.null(mask)) {
    mask = asNifti(mask) > 0
  } else {
    mask = array(TRUE, dim = dim(img))
  }

  if (packageVersion("neurobase") >= package_version("1.29.1")) {
    mat = neurobase::images_to_matrix(imgs, mask, verbose = verbose)
  } else {
    mat = neurobase::images2matrix(imgs, mask)
  }
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

