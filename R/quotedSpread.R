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
    #cl = min(8, parallel::detectCores()),
    windSpeed = windSpeedRas[],
    windDir = windDirRas[],
    meanDist = p[[1]],# rlnorm(1, log(p[[1]]), log(p[4])),#meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])), with estimating advectionDir
    maxDistance = maxDistance)
  return(out33)
})


#' @param massAttacksStack A rasterStack object with maps of number of attacked trees per pixel.
predictQuotedSpread <- function(massAttacksStack, # startYears, atksKnownNextYr,
                                dispersalKernel = "Weibull",
                                propPineRas, maxDistance, quotedSpread, p,
                                windDirStack, windSpeedStack, ...) {
  nams <- intersect(names(massAttacksStack), names(windDirStack))
  # nams <- names(massAttacksStack)
  names(nams) <- nams
  # lastOne <- length(nams)
  # nams <- nams# [-lastOne]
  massAttacksList <- raster::unstack(massAttacksStack[[nams]])
  windDirList <- raster::unstack(windDirStack[[nams]])
  windSpeedList <- raster::unstack(windSpeedStack[[nams]])
  propPineRas <- propPineRas
  reps <- 1

  # When we are in "predict mode", we don't know the attacks next year ...
  #   so we supply a map of pine
  whTos <- which(propPineRas[] >= 0.3)
  atksKnownNextYr <- data.table::data.table(pixels = whTos)

  out22 <- mcMap(nam = nams,
                 massAttacksRas = massAttacksList,
                 windSpeedRas = windSpeedList,
                 windDirRas = windDirList,
                 mc.cores = length(nams),
                 function(nam, massAttacksRas, windDirRas, windSpeedRas) {
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

