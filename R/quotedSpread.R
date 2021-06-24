quotedSpreadWeibull3 <- quote({
  to <- raster::xyFromCell(massAttacksRas, atksKnownNextYr$pixels)
  from <- raster::xyFromCell(massAttacksRas, starts)
  out33 <- SpaDES.tools::distanceFromEachPoint(
    to = to,
    from = from,
    angles = TRUE,
    cumulativeFn = "+",
    distFn = distanceFunction,
    landscape = massAttacksRas,
    propPineRas = propPineRas,
    dispersalKernel = dispersalKernel,
    p_advectionMag = p[2],
    p_sdDist = p[3],# p[4], # when there was estimation for p_advectionDir == needed p[4]
    cl = clNumber,
    windSpeed = windSpeedRas[],
    windDir = windDirRas[],
    p_meanDist = p[1],# rlnorm(1, log(p[[1]]), log(p[4])),#p_meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])), with estimating p_advectionDir
    maxDistance = maxDistance)
  return(out33)
})

quotedSpreadExponential <- quote({
  to <- raster::xyFromCell(massAttacksRas, atksKnownNextYr$pixels)
  from <- raster::xyFromCell(massAttacksRas, starts)
  out33 <- SpaDES.tools::distanceFromEachPoint(
    to = to,
    from = from,
    angles = TRUE,
    cumulativeFn = "+",
    distFn = distanceFunction,
    landscape = massAttacksRas,
    propPineRas = propPineRas,
    dispersalKernel = dispersalKernel,
    p_advectionMag = p[2],
    cl = clNumber,
    windSpeed = windSpeedRas[],
    windDir = windDirRas[],
    p_meanDist = p[1],# rlnorm(1, log(p[[1]]), log(p[4])),#p_meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])), with estimating p_advectionDir
    maxDistance = maxDistance)
  return(out33)
})

quotedSpreadGeneralGamma <- quote({
  to <- raster::xyFromCell(massAttacksRas, atksKnownNextYr$pixels)
  from <- raster::xyFromCell(massAttacksRas, starts)
  out33 <- SpaDES.tools::distanceFromEachPoint(
    to = to,
    from = from,
    angles = TRUE,
    cumulativeFn = "+",
    distFn = distanceFunction,
    landscape = massAttacksRas,
    propPineRas = propPineRas,
    dispersalKernel = dispersalKernel,
    p_advectionMag = p[2],
    p_sdDist = p[3],# p[4], # when there was estimation for p_advectionDir == needed p[4]
    cl = clNumber,
    windSpeed = windSpeedRas[],
    windDir = windDirRas[],
    p_meanDist = p[1],# rlnorm(1, log(p[[1]]), log(p[4])),#p_meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])), with estimating p_advectionDir
    p_nu = p[4],
    maxDistance = maxDistance)
  return(out33)
})


#' @param massAttacksStack A rasterStack object with maps of number of attacked trees per pixel.
predictQuotedSpread <- function(massAttacksDT,
                                dispersalKernel = "Weibull3", colNameForPrediction = "ATKTREES",
                                propPineRas, pineThreshold = 0.3, maxDistance, quotedSpread, p,
                                windDirStack, windSpeedStack, clNumber = NULL, ...) {
  if (!compareRaster(# massAttacksStack,
                     propPineRas,
                     windDirStack,
                     windSpeedStack))
    stop("The coarser resolution files need to be all the same resolution")

  nams <- intersect(unique(massAttacksDT$layerName), names(windDirStack))
  names(nams) <- nams
  windDirStk <- windDirStack[[nams]]
  windDirList <- if (nlayers(windDirStk) > 1) raster::unstack(windDirStk) else list(windDirStk)
  windSpeedStk <- windSpeedStack[[nams]]
  windSpeedList <- if (nlayers(windSpeedStk) > 1) raster::unstack(windSpeedStk) else list(windSpeedStk)
  propPineRas <- propPineRas

  # When we are in "predict mode", we don't know the attacks next year ...
  #   so we supply a map of pine
  whTos <- which(propPineRas[] >= pineThreshold)
  atksKnownNextYr <- data.table::data.table(pixels = whTos)

  out22 <- mcMap(nam = nams,
                 # massAttacksRas = massAttacksList,
                 windSpeedRas = windSpeedList,
                 windDirRas = windDirList,
                 MoreArgs = list(massAttacksDT),
                 mc.cores = length(nams),
                 function(nam, massAttacksDT, windDirRas, windSpeedRas) {
                   layerName1 <- names(windDirRas)
                   massAttacksDTThisYr <- massAttacksDT[layerName %in% layerName1]

                   massAttacksRas <- raster(windDirRas)
                   massAttacksRas[massAttacksDTThisYr$pixel] <- massAttacksDTThisYr$greenTreesYr_t # note green trees
                   starts <- which(massAttacksRas[] > 0)

                   if (length(starts) && length(whTos)) {
                     message("Predicting ", nam)
                     out <- as.data.table(eval(quotedSpread, envir = environment()))
                   } else {
                     out <- list()
                   }
                   return(out)
                 }
  )
  out22 <- rbindlist(out22, idcol = "layerName")
  out22[, `:=`(pixel = cellFromXY(propPineRas, cbind(x, y)))]
  out22[, predYear := yrNamesPlus1(layerName)]
  setnames(out22, "val", colNameForPrediction)
  out22
}

