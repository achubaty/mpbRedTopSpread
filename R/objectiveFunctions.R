
#' Objective function to minimize for mpbRedTopSpread
objFun <- function(p, propPineRas,
                   massAttacksStack, massAttacksDT,
                   quotedSpread,
                   omitPastPines, p_sdDist, dispersalKernel,
                   maxDistance, windDirStack, windSpeedStack) {

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
  if (missing(p_sdDist))
    p_sdDist <- get("p_sdDist", envir = .GlobalEnv)
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
    startYears = allYears[-length(allYears)],
    endYears = allYears[-1])

  maxObjFunValIndiv <- -Inf
  outBig <- purrr::pmap(
    .l = years, # this creates startYears and endYears because it is a data.frame
    massAttacksStack = massAttacksStack,
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
  startYears, endYears,
  massAttacksStack,
  massAttacksDT,
  clNumber = NULL,
  p,
  ...,
  # objFunValOnFail = 1e3,
  omitPastPines
) {
  layerName1 <- grep(startYears, names(massAttacksStack), value = TRUE)

  massAttacksDTThisYr <- massAttacksDT[layerName %in% layerName1]
  massAttacksRas <- raster(massAttacksStack[[startYears]])
  massAttacksRas[massAttacksDTThisYr$pixel] <- massAttacksDTThisYr$greenTreesYr_t # note green trees
  starts <- which(massAttacksRas[] > 0)
  windDirRas <- windDirStack[[startYears]]
  windSpeedRas <- windSpeedStack[[startYears]]
  atksRasNextYr <- massAttacksStack[[endYears]]
  # }
  wh <- which(atksRasNextYr[] > 0)
  atksKnownNextYr <- setDT(list(pixels = wh, ATKTREES = atksRasNextYr[][wh]))

  # this env has all the list(...), plus inherits the environment(), meaning all the objs & args above
  ll <- list2env(list(...))
  out <- # lapply(seq_len(ll$reps), function(rep)
    as.data.table(eval(quotedSpread, envir = ll))#)
  # out <- rbindlist(out, idcol = "rep")

  if (omitPastPines) {

    namesMAM <- names(massAttacksStack)

    namesToStartYears <- namesMAM[seq(which(namesMAM %in% startYears))]
    # unlist them
    mam <- if (length(namesToStartYears) == 1)
      list(massAttacksStack[[namesToStartYears]])
    else
      raster::unstack(massAttacksStack[[namesToStartYears]])
    names(mam) <- namesToStartYears
    massAttacksDTYearsToPres <- lapply(mam, function(x) {
      pixels = which(x[] > 0)
      setDT(list(pixels = pixels, ATKTREES = x[][pixels]))
    })
    massAttacksDTYearsToPres <- rbindlist(massAttacksDTYearsToPres, idcol = "layerName")
  }

  expectedNum <- out$val # * 1125 * prod(res(massAttacksRas))/1e4

  # Rescale -- we are interested in the distribution, not the absolute values -- maybe?: TODO
  # p_rescaler <- sum(atksKnownNextYr$ATKTREES)/sum(expectedNum)
  p_rescaler <- p["p_rescaler"]
  if (is.na(p_rescaler))
    p_rescaler <- 1 # 2.2 # was best in one case
  if (is.na(p_rescaler)) p_rescaler <- 1
  # likelihood
  lll <- dnbinom(round(atksKnownNextYr$ATKTREES), mu = expectedNum*p_rescaler, size = 0.2, log = TRUE)
  theInfs <- is.infinite(lll)
  if (any(theInfs)) {
    lowestProb <- min(lll[!theInfs])
    lll[theInfs] <- lowestProb
  }
  objFunVal <- lll

  objFunVal <- -round(sum(objFunVal, na.rm = TRUE), 3) # negative so optimizer does minimize
  # message(paste("Done startYear: ", startYears, " ObjFunVal: ", objFunVal))
  return(objFunVal)
}

