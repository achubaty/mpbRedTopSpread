
#' Objective function to minimize for mpbRedTopSpread
objFun <- function(p, propPineRas,
                   massAttacksStack, massAttacksDT,
                   quotedSpread,
                   omitPastPines, # p_sdDist,
                   dispersalKernel,
                   maxDistance, windDirStack, windSpeedStack,
                   objFunInner) {

  # DEoptim, when used across a network, inefficiently moves large objects every time
  #   it calls makes an objFun call. Better to move the objects to each node's disk,
  #   then read them in locally from disk --> faster when bandwidth is slow
  if (missing(propPineRas))
    propPineRas <- get("propPineRas", envir = .GlobalEnv)
  if (missing(massAttacksStack))
    massAttacksStack <- get("massAttacksStack", envir = .GlobalEnv)
  if (missing(massAttacksDT))
    massAttacksDT <- get("massAttacksDT", envir = .GlobalEnv)
  if (missing(omitPastPines))
    omitPastPines <- get("omitPastPines", envir = .GlobalEnv)
  # if (missing(p_sdDist))
  #   p_sdDist <- get("p_sdDist", envir = .GlobalEnv)
  if (missing(dispersalKernel))
    dispersalKernel <- get("dispersalKernel", envir = .GlobalEnv)
  if (missing(maxDistance))
    maxDistance <- get("maxDistance", envir = .GlobalEnv)
  if (missing(windDirStack))
    windDirStack <- get("windDirStack", envir = .GlobalEnv)
  if (missing(windSpeedStack))
    windSpeedStack <- get("windSpeedStack", envir = .GlobalEnv)

  allYears <- names(massAttacksStack)
  years <- data.frame(
    t0 = allYears[-length(allYears)],
    tEnd = allYears[-1],
    tminus1 = names(windDirStack))

  # sampleSizes <- massAttacksDT[, .N, by = "layerName"]
  # keepLayers <- sampleSizes$layerName[sampleSizes$N > 1]
  # years <- years[years$t0 %in% keepLayers & years$tEnd %in% keepLayers,]
  # years$t0
  # years <- years[years$t0 != "X2004", ] # has 1 origin pixel
  maxObjFunValIndiv <- -Inf

  outBig <- purrr::pmap(
    .l = years, # this creates t0 and tEnd because it is a data.frame
    rasterTemplate = terra::rast(massAttacksStack[[1]]),
    # massAttacksStack = massAttacksStack,
    massAttacksDT = massAttacksDT,
    maxDistance = maxDistance,
    p = p, # reps = reps,
    propPineRas = propPineRas,
    quotedSpread = quotedSpread,
    omitPastPines = omitPastPines,
    windDirStack = windDirStack,
    windSpeedStack = windSpeedStack,
    dispersalKernel = dispersalKernel,
    objFunInner # the function
  )

  sum(unlist(outBig), na.rm = TRUE)

}

#' @param ... Any arguments that are used only in the quotedSpread
objFunInner <- function(
  quotedSpread,
  t0, tEnd, tminus1,
  rasterTemplate,
  # massAttacksStack,
  massAttacksDT,
  clNumber = NULL,
  p, toCells,
  ...,
  # objFunValOnFail = 1e3,
  omitPastPines,
  mode = c("fit", "predict")
) {

  mode <- mode[1]
  layerNames <- unique(massAttacksDT$layerName)
  t0LayerNam <- grep(t0, layerNames, value = TRUE)
  tminus1LayerNam <- grep(tminus1, names(windDirStack), value = TRUE)


  massAttacksDTThisYr <- massAttacksDT[layerName %in% t0LayerNam]
  # rasterTemplate <- rast(massAttacksStack[[t0]])
  massAttackRasStart <- rasterTemplate
  massAttackRasStart[massAttacksDTThisYr$pixel] <- massAttacksDTThisYr$greenTreesYr_t # note green trees
  totalPropaguleNumber <- sum(massAttacksDTThisYr$greenTreesYr_t)
  starts <- which(massAttackRasStart[] > 0)
  windDirRas <- windDirStack[[tminus1LayerNam]]
  windSpeedRas <- windSpeedStack[[tminus1LayerNam]]

  # tEnd
  madtTend <- massAttacksDT[layerName %in% tEnd]# $ATKTREES
  if (missing(toCells)) {
    toCells <- whPixelsNextYear
  }
  ATKTREESNextYear <- madtTend$ATKTREES
  # atksRasNextYr <- massAttacksStack[[tEnd]]
  # }
  # wh <- which(atksRasNextYr[] > 0)
  # atksKnownNextYr <- setDT(list(pixels = whPixelsNextYear, ATKTREES = ATKTREESNextYear))

  # this env has all the list(...), plus inherits the environment(), meaning all the objs & args above

  ll <- list2env(list(...))
  out <- as.data.table(eval(quotedSpread, envir = ll))#)

  if (omitPastPines) {

    # layerNames <- names(massAttacksStack)

    namesToStartYears <- layerNames[seq(which(layerNames %in% t0))]
    # unlist them
    mam <- if (length(namesToStartYears) == 1)
      list(massAttacksStack[[namesToStartYears]])
    else
      as.list(massAttacksStack[[namesToStartYears]])
    names(mam) <- namesToStartYears
    massAttacksDTYearsToPres <- lapply(mam, function(x) {
      pixels = which(x[] > 0)
      setDT(list(pixels = pixels, ATKTREES = x[][pixels]))
    })
    massAttacksDTYearsToPres <- rbindlist(massAttacksDTYearsToPres, idcol = "layerName")
  }

  browser()
  expectedNum <- out$val # * 1125 * prod(res(massAttackRasStart))/1e4

  # Rescale -- we are interested in the distribution, not the absolute values -- maybe?: TODO
  p_rescaler <- if (mode == "fit") p_rescaler <- sum(ATKTREESNextYear)/totalPropaguleNumber else 18

  # p_rescaler <- sum(ATKTREESNextYear)/sum(expectedNum)
  # p_rescaler <- p["p_rescaler"]
  # if (is.na(p_rescaler))
  #   p_rescaler <- 1 # 2.2 # was best in one case
  # if (is.na(p_rescaler)) p_rescaler <- 1

  # likelihood
  if (identical(mode, "fit")) {
    lll <- dnbinom(round(ATKTREESNextYear), mu = expectedNum*p_rescaler, size = 0.2, log = TRUE)
    theInfs <- is.infinite(lll)
    if (any(theInfs)) {
      if (all(theInfs)) {
        lowestProb <- -Inf
      } else {
        lowestProb <- min(lll[!theInfs])
      }
      lll[theInfs] <- lowestProb
    }
    objFunVal <- lll

    out <- -round(sum(objFunVal, na.rm = TRUE), 3) # negative so optimizer does minimize
  } else {
    out <- expectedNum*p_rescaler
    # out <- rnbinom(n = length(expectedNum), mu = expectedNum*p_rescaler, size = 0.2)
  }
  return(out)
}
