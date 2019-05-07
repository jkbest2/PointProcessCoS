
## Function to generate a Polygon of a block from the lower-left corner
block_poly <- function(ll, w = 25) {
  coords <- matrix(c(ll[1],     ll[2],
                     ll[1] + w, ll[2],
                     ll[1] + w, ll[2] + w,
                     ll[1],     ll[2] + w,
                     ll[1],     ll[2]),
                   byrow = TRUE, ncol = 2L)
  Polygon(coords, hole = FALSE)
}

## Get elements of intensity that are within a block. Takes a Polygon and the
## saved "Lambda" attribute of the log-Gaussian Cox process generated above.
## This is a rough integral over each block. The `eps` argument in the `rLGCP`
## call above controls the grid spacing and so the integral approximation.
total_intens <- function(blk_poly, Lambda) {
  crds <- blk_poly@coords
  x_idx <- which((Lambda$xcol > crds[1, 1]) & (Lambda$xcol < crds[3, 1]))
  y_idx <- which((Lambda$yrow > crds[1, 2]) & (Lambda$yrow < crds[3, 2]))
  sum(Lambda$v[x_idx, y_idx])
}

## Calculate the weights (areas) of each polygon of the dual mesh using `rgeos`
## functions, intersected with the point process domain. Taken from Chapter 4 of
## Advanced Spatial Modeling (linked above).
get_wts <- function(poly, mesh_dual) {
  wts <- vapply(1:length(mesh_dual), function(idx) {
    if (gIntersects(mesh_dual[idx, ], poly)) {
      area <- gArea(gIntersection(mesh_dual[idx, ], poly))
    } else {
      area <- 0
    }
    area
  }, 1.0)
  wts
}

count_points <- function(points, mesh_dual, mask = pp_dom_poly) {
  ## Use mask to count only points withing mask polygons. In this case we use
  ## the chosen quadrats. Otherwise count within the whole domain.
  vapply(1:length(mesh_dual), function(idx) {
    ## This intersects the dual mesh with the mask (quadrats), the intersects
    ## that small polygon with the point process. Finally, it counts the number
    ## of points that fall within that polygon, which is what it returns.
    dual_poly <- gIntersection(mesh_dual[idx, ], mask)
    if (is.null(dual_poly)) {
      0L
    } else {
      length(gIntersection(dual_poly, pp_sp))
    }
  }, 1L)
}
