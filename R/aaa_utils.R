remake_nifti_image = function(vec, img, mask) {
  vec = c(vec)
  arr = array(NA_real_, dim = dim(mask))
  arr[ mask > 0 ] = vec
  arr = RNifti::asNifti(arr, reference = img)
  arr
}


