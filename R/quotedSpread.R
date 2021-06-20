quotedSpreadDefault <- quote({
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
    asymParam = p[2],
    sdDist = p[3],# p[4], # when there was estimation for advectionDir == needed p[4]
    cl = clNumber,
    windSpeed = windSpeedRas[],
    windDir = windDirRas[],
    meanDist = p[[1]],# rlnorm(1, log(p[[1]]), log(p[4])),#meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])), with estimating advectionDir
    maxDistance = maxDistance)
  return(out33)
})


#' @param massAttacksStack A rasterStack object with maps of number of attacked trees per pixel.
predictQuotedSpread <- function(massAttacksDT, massAttacksStack, # startYears, atksKnownNextYr,
                                dispersalKernel = "Weibull",
                                propPineRas, pineThreshold = 0.3, maxDistance, quotedSpread, p,
                                windDirStack, windSpeedStack, clNumber = NULL, ...) {
  if (!compareRaster(massAttacksStack,
                     propPineRas,
                     windDirStack,
                     windSpeedStack))
    stop("The coarser resolution files need to be all the same resolution")

  nams <- intersect(unique(massAttacksDT$layerName), names(windDirStack))
  # nams <- names(massAttacksStack)
  names(nams) <- nams
  # lastOne <- length(nams)
  # nams <- nams# [-lastOne]
  # massAttacksList <- raster::unstack(massAttacksStack[[nams]])
  windDirStk <- windDirStack[[nams]]
  windDirList <- if (nlayers(windDirStk) > 1) raster::unstack(windDirStk) else list(windDirStk)
  windSpeedStk <- windSpeedStack[[nams]]
  windSpeedList <- if (nlayers(windSpeedStk) > 1) raster::unstack(windSpeedStk) else list(windSpeedStk)
  propPineRas <- propPineRas
  reps <- 1

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
                     out <- lapply(seq_len(reps), function(rep)
                       as.data.table(eval(quotedSpread, envir = environment())))
                     out <- rbindlist(out, idcol = "rep")
                   } else {
                     out <- list()
                   }
                   return(out)
                 }
  )
  rbindlist(out22, idcol = "Year")
}

