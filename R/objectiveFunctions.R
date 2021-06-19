
objFun <- function(p, propPineRas, # massAttacksRas,
                   advectionDir, advectionMag,
                   # minNumAgents,
                   massAttacksStack, massAttacksDT, currentTime, reps = 10,
                   quotedSpread, fitType = "ss1", omitPastPines, sdDist, dispersalKernel,
                   maxDistance, windDirStack, windSpeedStack, subsampleFrom = NULL) {

  # DEoptim, when used across a network, inefficiently moves large objects every time
  #   it calls makes an objFun call. Better to move the objects to each node's disk,
  #   then read them in locally from disk --> faster when bandwidth is slow
  if (missing(propPineRas))
    propPineRas <- get("propPineRas", envir = .GlobalEnv)
  if (missing(massAttacksStack))
    massAttacksStack <- get("massAttacksStack", envir = .GlobalEnv)
  if (missing(massAttacksDT))
    massAttacksDT <- get("massAttacksDT", envir = .GlobalEnv)
  # if (missing(massAttacksRas))
  #   massAttacksRas <- get("massAttacksRas", envir = .GlobalEnv)
  if (missing(advectionDir))
    advectionDir <- get("advectionDir", envir = .GlobalEnv)
  if (missing(advectionMag))
    advectionMag <- get("advectionMag", envir = .GlobalEnv)
  # if (missing(minNumAgents))
  #   minNumAgents <- get("minNumAgents", envir = .GlobalEnv)
  if (missing(currentTime))
    currentTime <- get("currentTime", envir = .GlobalEnv)
  if (missing(fitType))
    fitType <- get("fitType", envir = .GlobalEnv)
  if (missing(omitPastPines))
    omitPastPines <- get("omitPastPines", envir = .GlobalEnv)
  if (missing(sdDist))
    sdDist <- get("sdDist", envir = .GlobalEnv)
  if (missing(dispersalKernel))
    dispersalKernel <- get("dispersalKernel", envir = .GlobalEnv)
  if (missing(maxDistance))
    maxDistance <- get("maxDistance", envir = .GlobalEnv)
  if (missing(windDirStack))
    windDirStack <- get("windDirStack", envir = .GlobalEnv)
  if (missing(windSpeedStack))
    windSpeedStack <- get("windSpeedStack", envir = .GlobalEnv)
  if (missing(subsampleFrom))
    subsampleFrom <- get("subsampleFrom", envir = .GlobalEnv)

  if (fitType == "likelihood") {
    if (reps < 10) {
      message("reps must be at least 10 if using likelihood; setting to 10")
      reps <- 10
    }
  }

  allYears <- names(massAttacksStack)

  combinations <- seq_along(allYears[-1])
  years <- data.frame(
    startYears = allYears[-length(allYears)],
    endYears = allYears[-1])
  # combinations <- cbind(combinations, years[combinations[, "years"], ])
  # combinations <- as.data.frame(combinations[, c("reps", "startYears", "endYears")])

  maxObjFunValIndiv <- -Inf
  # mam <- raster::unstack(massAttacksStack)
  outBig <- purrr::pmap(
    .l = years, # this creates startYears and endYears because it is a data.frame
    massAttacksStack = massAttacksStack,
    massAttacksDT = massAttacksDT,
    maxDistance = maxDistance,
    p = p, reps = reps,
    # minNumAgents = minNumAgents,
    propPineRas = propPineRas,
    quotedSpread = quotedSpread,
    fitType = fitType,
    omitPastPines = omitPastPines,
    windDirStack = windDirStack,
    windSpeedStack = windSpeedStack,
    dispersalKernel = dispersalKernel,
    subsampleFrom = subsampleFrom,
    .f = objFunInner
  )

  sum(unlist(outBig), na.rm = TRUE)

}

#' @param ... Any arguments that are used only in the quotedSpread
objFunInner <- function(
  quotedSpread,
  startYears, endYears,
  massAttacksStack,
  massAttacksDT,
  ...,
  objFunValOnFail = 1e3,
  fitType, omitPastPines,
  subsampleFrom = NULL
) {
  layerName1 <- grep(startYears, names(massAttacksStack), value = TRUE)

  massAttacksDTThisYr <- massAttacksDT[layerName %in% layerName1]
  if (!is.null(subsampleFrom)) {
    if (Require:::isWindows() || amc::isRstudio())
      message("subsampleFrom is not NULL; using ", floor(NROW(massAttacksDTThisYr)/subsampleFrom),
              " from cells, about 1/", subsampleFrom)
    massAttacksDTThisYr <- massAttacksDTThisYr[sort(sample(1:NROW(massAttacksDTThisYr), NROW(massAttacksDTThisYr)/subsampleFrom))]
  }
  # starts <- massAttackThisYr$pixel
  # windDirVec <- massAttackThisYr$windDir
  # windSpeedVec <- massAttackThisYr$windSpeed
  # massAttackNextYr <- massAttackThisYr$ATKTREEStplus1
  # if (FALSE) {
  # massAttacksRas <- massAttacksStack[[startYears]]
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
  out <- lapply(seq_len(ll$reps), function(rep) as.data.table(eval(quotedSpread, envir = ll)))
  out <- rbindlist(out, idcol = "rep")

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
    massAttacksDTYearsToPres <- rbindlist(massAttacksDTYearsToPres, idcol = "Year")
  }
  if (fitType == "distanceFromEachPoint") {
    expectedNum <- out$val # * 1125 * prod(res(massAttacksRas))/1e4

    # Rescale -- we are interested in the distribution, not the absolute values -- maybe?: TODO
    # rescaler <- sum(atksKnownNextYr$ATKTREES)/sum(expectedNum)
    rescaler <- 1#sum(atksKnownNextYr$ATKTREES)/sum(expectedNum)
    # likelihood
    lll <- dnbinom(round(atksKnownNextYr$ATKTREES), mu = expectedNum*rescaler, size = 0.1, log = TRUE)
    theInfs <- is.infinite(lll)
    if (any(theInfs)) {
      lowestProb <- min(lll[!theInfs])
      lll[theInfs] <- lowestProb
    }
    objFunVal <- lll

  } else {
    stop("not tested recently; use fitType = 'distanceFromEachPoint'")
    # Remove pixels that had already been attacked in the past -- emulating MPB Suppression efforts
    if (omitPastPines)
      atksKnownNextYr <- atksKnownNextYr[!massAttacksDTYearsToPres, on = "pixels"]

    atksNextYearSims <- out[, list(abundSettled  = sum(abundSettled) ), by = c("rep", "pixels")]
    # nPix <- atksNextYearSims[abundSettled > 0, .N] ## total number of pixels

    ## attacked area from data
    # atksKnownNextYr <- atksKnownNextYr[ATKTREES > 0]
    # i <- 1
    oo <- atksKnownNextYr[atksNextYearSims, on = "pixels", nomatch = NA]
    set(oo, which(is.na(oo$ATKTREES)), "ATKTREES", 0L)
    # oo <- atksNextYearSims[atksKnownNextYr, on = "pixels", nomatch = NA]
    if (reps > 1)
      oo[, ATKTREES := ATKTREES[1], by = "pixels"]

    #i <- 0
    if (fitType == "likelihood") {
      stop("not tested recently")
      probs <- oo[, list(prob = {
        prob <- if (.N > 2) {
          i <<- i + 1
          # print(i)

          max(1e-14, demp(ATKTREES[1], abundSettled))

        } else {
          1e-14
        }
      }), by = "pixels"]
      objFunVal <- sum(-log(probs$prob))
    } else {
      oo[is.na(abundSettled),  abundSettled := 0]
      if (fitType == "ss1") {
        objFunVal <- sum((oo$ATKTREES - oo$abundSettled)^2)
      } else {
        #if (startYears == "X2010") browser()
        # To balance the zeros and non-zeros, must sub-sample the zeros so there are equal
        #   number as the non-zeros
        whNonZero <- oo$ATKTREES > 0
        numNonZero <- tabulate(as.integer(whNonZero))
        whZero <- which(!whNonZero)
        whNonZero <- which(whNonZero)

        # There are way more zeros than non zeros. We need to sub-sample the zeros or
        #   they influence the objFun too much. Here, hard coded -- take 1/4 of the number
        #   of zeros as there are non-zeros
        samZero <- sample(whZero, size = round(numNonZero/4, 0))
        samNonZero <- sample(whNonZero, size = numNonZero)
        sam <- c(samZero, samNonZero)
        atkTrees <- oo$ATKTREES[sam]
        abundSett <- oo$abundSettled[sam]

        # Rescale so that this is a relative abundance
        ratio <- sum(abundSett)/sum(atkTrees)
        abundSettRescaled <- abundSett/ratio

        if (fitType == "logSAD") {
          objFunVal <- log(abs(atkTrees - abundSettRescaled))
        } else if (fitType == "SAD") {
          objFunVal <- abs(atkTrees - abundSettRescaled)
        } else {
          stop("fitType must be logSAD or SAD")
        }
        # browser()
        isInf <- is.infinite(objFunVal)
        if (any(isInf)) {
          if (all(isInf)) {
            negInf <- objFunVal < 0
            if (any(negInf))
              objFunVal[isInf] <- 0
            if (any(!negInf))
              objFunVal[isInf] <- objFunValOnFail
          } else {
            objFunVal[isInf] <- max(objFunVal[!isInf], na.rm = TRUE)
          }
        }
      }
    }
  }
  objFunVal <- -round(sum(objFunVal, na.rm = TRUE), 3) # negative so optimizer does minimize
  # message(paste("Done startYear: ", startYears, " ObjFunVal: ", objFunVal))
  return(objFunVal)
}

