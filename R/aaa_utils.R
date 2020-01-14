remake_nifti_image = function(vec, img, mask) {
  vec = c(vec)
  arr = array(NA_real_, dim = dim(mask))
  arr[ mask ] = vec
  arr = RNifti::asNifti(image = arr, reference = img)
  arr
}


