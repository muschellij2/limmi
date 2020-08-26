#' Compute local moderated T-statistics
#'
#' @param imgs List of images or character of filenames
#' @param mask Image to subset the data
#' @param verbose Print diagnostic messages
#' @param nn nearest-neighbors computed by \code{\link{findNearestNeighbors}}
#' @param radius Integer: Indices that dist(targetInd - queryInd) <= radius
#' @param  adjust.method method to adjust p-value,
#'  passed to \code{\link{topTable}}
#' @param ... additional arguments to pass to \code{\link[parallel]{mclapply}}

#' @return A list of output images, including the moderated T map, p-value map,
#' adjusted p-value map, and the pseduo T map.
#' @export


nifti_local_moderated_t_stat = function(imgs, mask = NULL, nn, radius, adjust.method = "BH", verbose = FALSE, ...){

  maps = nifti_images_to_matrix(imgs = imgs, mask = mask)
  nrow.map = nrow(maps)

  t.stat.list = parallel::mclapply(1:nrow.map,
                     function(center){
                       # require(limma)
                       idx = nn$index[center, nn$distance[center,]<= radius]
                       fit = limma::lmFit(maps[idx,])
                       eb.fit = limma::eBayes(fit)
                       tstat = eb.fit$t[1]
                       pval = eb.fit$p.value[1]
                       # local.eb.fit = sapply(eb.fit, "[", 1)
                       adj.pval = limma::topTable(eb.fit, sort.by = 'none', adjust.method = adjust.method, number = length(idx))$adj.P.Val[1]
                       psuedo.t = eb.fit$coefficients[1]/sqrt(mean((eb.fit$sigma)^2))/eb.fit$stdev.unscaled[1]
                       return(c(tstat, pval, adj.pval, psuedo.t))
                     })

  t.stat = do.call("cbind", t.stat.list)

  tmap.limma = mask * NA
  pvalmap.limma = mask * NA
  adj.pvalmap.limma = mask * NA
  pseudo.tmap.limma = mask * NA
  if (is.null(mask)){
    tmap.limma = t.stat[1,]
    pvalmap.limma = t.stat[2,]
    adj.pvalmap.limma = t.stat[3,]
    pseudo.tmap.limma = t.stat[4,]
  }else{
    tmap.limma[mask>0] = t.stat[1,]
    pvalmap.limma[mask>0] = t.stat[2,]
    adj.pvalmap.limma[mask>0] = t.stat[3,]
    pseudo.tmap.limma[mask>0] = t.stat[4,]
  }
  return(list(tmap = tmap.limma,
              pvalmap = pvalmap.limma,
              adj.pvalmap = adj.pvalmap.limma,
              pseudo.tmap = pseudo.tmap.limma)
         )

}

