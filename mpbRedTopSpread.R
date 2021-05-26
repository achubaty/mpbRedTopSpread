defineModule(sim, list(
  name = "mpbRedTopSpread",
  description = "Mountain Pine Beetle Red Top Growth Model: Short-run Potential for Establishment, Eruption, and Spread",
  keywords = c("mountain pine beetle, outbreak dynamics, eruptive potential, spread, climate change, twitch response"),
  authors = c(
    person(c("Alex", "M"), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre")),
    person(c("Eliot", "J B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut")),
    person(c("Barry", "J"), "Cooke", email = "barry.cooke@canada.ca", role = c("aut"))
  ),
  childModules = character(),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  reqdPkgs = list("achubaty/amc@development", "CircStats", "data.table",
                  "DEoptim", "EnvStats",
                  "PredictiveEcology/LandR@development", "parallelly",
                  "PredictiveEcology/pemisc@development (>= 0.0.3.9001)",
                  "quickPlot", "raster", "RColorBrewer", "reproducible",
                  "PredictiveEcology/SpaDES.tools@spread3 (>= 0.3.7.9018)"),
  parameters = rbind(
    defineParameter("advectionDir", "numeric", 90, 0, 359.9999,
                    "The direction of the spread bias, in degrees from north"),
    defineParameter("advectionMag", "numeric", 1000, NA, NA,
                    "The magnitude of the directional bias of spread"),
    defineParameter("bgSettlingProp", "numeric", 0.1, 0, 1,
                    "The proportion of beetles that settle from those that could potentially settle, even if no pine"),
    defineParameter("meanDist", "numeric", 1000, NA, NA,
                    "Expected dispersal distance (m); ~63% go less than this distance"),
    defineParameter("type", "character", "fit", NA, NA,
                    "Expected dispersal distance (m); ~63% go less than this distance"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "This describes the interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should this entire module be run with caching activated?")
  ),
  inputObjects = bindrows(
    expectsInput("currentAttacks", "RasterLayer",
                 desc = "Current MPB attack map (number of red attacked trees).",
                 sourceURL = NA),
    expectsInput("massAttacksDT", "data.table",
                 desc = "Current MPB attack map (number of red attacked trees).",
                 sourceURL = NA),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "if not supplied, will default to standAgeMap", # TODO: description needed
                 sourceURL = NA),
    expectsInput("standAgeMap", "RasterLayer",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
  ),
  outputObjects = bindrows(
    createsOutput("massAttacksDT", "data.table", "Current MPB attack map (number of red attacked trees).")
  )
))

## event types
#   - type `init` is required for initilization

doEvent.mpbRedTopSpread <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
         "init" = {
           # do stuff for this event
           sim <- Init(sim)

           # schedule future event(s)
           sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "dispersal", eventPriority = 4.5)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "plot", .last() - 1)
           sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbRedTopSpread", "save", .last())
         },
         "dispersal" = {
           out <- dispersal2(pineMap = sim$pineMap, studyArea = sim$studyArea,
                             massAttacksDT = sim$massAttacksDT,
                             massAttacksMap = sim$massAttacksMap,
                             currentAttacks = sim$currentAttacks,
                             params = P(sim),
                             bgSettlingProp = P(sim)$bgSettlingProp,
                             type = P(sim)$type,
                             currentTime = time(sim),
                             reqdPkgs = reqdPkgs(module = currentModule(sim), modulePath = modulePath(sim))[[currentModule(sim)]]
           )
           # sim <- dispersal(sim)
           sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "dispersal", eventPriority = 4.5)
         },
         "plot" = {
           sim <- plotFn(sim)
           sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "plot", eventPriority = .last() - 1)
         },
         warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                       "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  ## raster to match
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- Cache(
      LandR::prepInputsLCC,
      year = 2005,
      destinationPath = dPath,
      studyArea = sim$studyArea
    )
  }

  ## stand age map
  if (!suppliedElsewhere("standAgeMap", sim)) {
    sim$standAgeMap <- LandR::prepInputsStandAgeMap(
      startTime = 2010,
      ageUrl = na.omit(extractURL("standAgeMap")),
      destinationPath = dPath,
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      userTags = c("stable", currentModule(sim)) ## TODO: does this need rasterToMatch? it IS rtm!
    )
    sim$standAgeMap[] <- asInteger(sim$standAgeMap[])
  }

  return(invisible(sim))
}

### initilization
Init <- function(sim) {
  cr <- compareRaster(sim$rasterToMatch, sim$currentAttacks, sim$pineMap, orig = TRUE)
  if (!isTRUE(cr))
    stop("Rasters sim$rasterToMatch, sim$currentAttacks, sim$pineMap must all have same metadata; they do not")

  return(invisible(sim))
}

### plotting
plotFn <- function(sim) {
  r <- raster(sim$massAttacksMap[[1]]) ## use a template
  atkMap <- amc::dt2raster(sim$massAttacksDT, r, "ATKTREES")
  setColors(atkMap) <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
  Plot(atkMap, title = "Simulated Attacks")

  return(invisible(sim))
}

### spread
dispersal2 <- function(pineMap, studyArea, massAttacksDT, massAttacksMap,
                       currentAttacks, params, currentTime, bgSettlingProp, type, reqdPkgs) {
  ## check that MPB and pine rasters are the same resolution and ncells
  if (fromDisk(pineMap))
    pineMap[] <- pineMap[]

  ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (2017)
  MAXTREES <- 1125 * prod(res(pineMap)) / 100^2 ## TODO: round this?

  ## asymmetric spread (biased eastward)
  # lodgepole pine and jack pine together
  propPineMap <- if (nlayers(pineMap) > 1) sum(pineMap) else pineMap
  mv <- maxValue(propPineMap)

  if (mv > 1)
    propPineMap[] <- propPineMap[] / 100 ## much faster than calc; ## TODO: allow different params per species

  nas <- is.na(propPineMap[])
  propPineMap[nas] <- 0
  if (exists("EliotTesting")) {
    EliotTesting <- TRUE
    sim@params$mpbRedTopSpread$bgSettlingProp <- 0.7
    sim@params$mpbRedTopSpread$advectionMag <- 20000
    minNumAgents <- 200
    sim@params$mpbRedTopSpread$.plotInitialTime <- NA
  } else {
    EliotTesting <- FALSE
  }
  propPineMap[] <- pmin(1, propPineMap[] + bgSettlingProp) ## TODO: why divide by 100??

  if (EliotTesting) { # TODO -- delete EliotTesting when no longer desired
    a <- extent(studyArea)
    starts <- massAttacksDT[ATKTREES > 0]$ID
    d <- raster(currentAttacks)
    d[starts] <- 1
    d <- crop(d, a)
    starts <- which(d[] > 0)
    currentAttacks <- crop(currentAttacks, a) * 300
    propPineMap <- crop(propPineMap, a)
    saveStack <- raster::rasterTmpFile()
  } else {
    minNumAgents <- 50
    starts <- massAttacksDT[["ID"]][massAttacksDT$ATKTREES > 0]
    saveStack <- NULL
    currentAttacks <- currentAttacks
  }

  # mam <- raster::unstack(massAttacksMap)
  # names(mam) <- names(massAttacksMap)
  # massAttacksDTAllYears <- lapply(mam, function(x) {
  #   pixels = which(x[] > 0)
  #   setDT(list(pixels = pixels, ATKTREES = x[][pixels]))
  # })

  # Put objects in Global for objFun -- this is only when not using multi-machine cluster
  advectionDir <- params$advectionDir
  advectionMag <- params$advectionMag
  omitPastPines <- TRUE
  sdDist <- 1.2
  dispersalKernel <- "Weibull"
  # dispersalKernel <- "Exponential"
  p <- do.call(c, params[c("meanDist", "advectionMag", "advectionDir")])
  p["meanDist"] <- 1e4
  p["advectionMag"] <- 5e3
  p["advectionDir"] <- 90
  p["sdDist"] <- 1.2
  p["meanDistSD"] <- 20
  p["advectionDirSD"] <- 20
  lower <- c(25000,  3000,   0, 0.9, 1.001,  5)
  upper <- c(45000, 20000, 330, 2.5, 1.6, 30)
  p[] <- sapply(seq_along(p), function(x) runif(1, lower[x], upper[x]))

  fitType <- "distanceFromEachPoint"
  objsToExport <- setdiff(formalArgs("objFun"), c("p", "reps", "quotedSpread", "fitType", "distanceFunction"))
  libPaths <- .libPaths()
  objsToExport <- c("reqdPkgs", objsToExport, "objsToExport", "libPaths")
  list2env(mget(objsToExport), envir = .GlobalEnv)

  quotedSpread <-
    quote(SpaDES.tools::spread3(start = starts,
                                rasQuality = propPineMapInner,
                                rasAbundance = currentAttacks,
                                advectionDir = rnorm(1, p[3], p[6]),
                                advectionMag = p[2],
                                meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])),
                                sdDist = p[4],
                                dispersalKernel = dispersalKernel,
                                plot.it = FALSE,
                                minNumAgents = minNumAgents,
                                verbose = 0,
                                skipChecks = TRUE,
                                saveStack = NULL))


  quotedSpread <-
    quote(distanceFromEachPoint(
      to = xyFromCell(atksRasNextYr, atksKnownNextYr$pixels),
      from = xyFromCell(atksRasNextYr, starts),
      angles = T, cumulativeFn = "+",
      distFn = distanceFunction,
      landscape = currentAttacks,
      nextYrVec = atksRasNextYr[],
      propPineMapInner = propPineMapInner,
      dispersalKernel = dispersalKernel,
      asym = p[2],
      sdDist = p[4],
      asymDir = rnorm(1, p[3], p[6]),
      meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])),
      maxDistance = 7e4))

  browser()
  if (isTRUE(type == "fit")) {
     # cl <- makeOptimalCluster(type = "FORK", MBper = 3000, min(3, length(p) * 10),
     #                          assumeHyperThreads = TRUE)

    # browser()
    message("Starting cluster with 1 core per machine -- install packages; copy objects; write to disk")
    clusterIPs <- c(rep("localhost", 26), rep("10.20.0.184", 12), rep("10.20.0.97", 22))
    clSingle <- future::makeClusterPSOCK(workers = unique(clusterIPs), revtunnel = TRUE)
    on.exit(try(parallel::stopCluster(clSingle)))
    clusterExport(clSingle, varlist = objsToExport, envir = environment())
    clusterEvalQ(clSingle, {
      if (any(!dir.exists(libPaths)))
        dir.create(libPaths[1], recursive = TRUE)
      .libPaths(libPaths)
      ## set CRAN repos; use binary linux packages if on Ubuntu
      local({
        options(Ncpus = parallel::detectCores() / 2)
        options("repos" = c(CRAN = "https://cran.rstudio.com"))

        if (Sys.info()["sysname"] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
          .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
          .user.agent <- paste0(
            "R/", getRversion(), " R (",
            paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"]),
            ")"
          )
          options(repos = c(CRAN = paste0("https://packagemanager.rstudio.com/all/__linux__/",
                                          .os.version, "/latest")))
          options(HTTPUserAgent = .user.agent)
        }
      })

      if (!suppressWarnings(require("Require"))) {
        install.packages("Require")
      }

      suppressMessages(Require::Require(reqdPkgs, install = TRUE, require = FALSE))
      save(list = objsToExport, file = "allObjs.rda")
    })
    parallel::stopCluster(clSingle)
    message("Starting cluster with all cores per machine")
    cl <- future::makeClusterPSOCK(workers = clusterIPs, revtunnel = TRUE)
    on.exit(parallel::stopCluster(cl))
    st <- system.time(clusterEvalQ(cl, {
      load(file = "allObjs.rda", envir = .GlobalEnv)
      .libPaths(libPaths)
      suppressMessages(Require::Require(c("SpaDES.tools", "raster", "CircStats", "data.table"), install = FALSE))
    }))
    message("Starting DEoptim")
    stPre <- Sys.time()
    DEout <- DEoptim(fn = objFun,
                     lower = lower,#c(500, 1, -90, 0.9, 1.1, 5),
                     upper = upper,#c(30000, 10, 240, 1.8, 1.6, 30),
                     reps = 1,
                     quotedSpread = quotedSpread,
                     control = DEoptim.control(cluster = cl),
                     fitType = fitType)
    saveRDS(DEout, file = file.path("outputs", paste0("DEout_", format(stPre), ".rds")))
    stPost <- Sys.time()
    print(difftime(stPost, stPre))

  } else {
    profvis::profvis(
      out <- objFun(quotedSpread = quotedSpread, reps = 1, p = p, fitType = fitType,
                  massAttacksMap = massAttacksMap)
    )
  }
  browser()
  objs <- dir("outputs", pattern = "DEout", full.names = TRUE)
  DEout <- readRDS(tail(objs, 1))
  # p <- numeric()
  # p[] <- c(38246.790564,   41.199916,  237.629198,    2.297772,    1.042395,   18.800792)
  p[] <- c(28386.900711,15887.312435,  226.518642,    1.602424 ,   1.557437,   26.188510)
  # p[] <- DEout$optim$bestmem[] # use [] to keep names


  N1 <- 1e4

  par(mfrow = c(3,4))
  for (i in 1:10) {
    meanDist = rlnorm(1, log(p[[1]]), log(p[[5]]))
    sdDist = p[4]

    mn <- (meanDist)
    sd <- mn/sdDist # 0.8 to 2.0 range
    shape <- (sd/mn)^(-1.086)
    scale <- mn/exp(lgamma(1+1/shape))

    ww <- rweibull(1e4, shape = shape, scale = scale)
    hist(ww, xlim = c(0, 14e4), main = "", xlab = "Distance (m) from source")
    abline(v = meanDist, col = "red")
    abline(v = p[[1]], col = "green")

  }
  maxD <- 70000 + 125
  r <- raster(extent(-maxD+ 250, maxD + 250 , -maxD, maxD), res = 250)
  r[] <- 1
  r <- gaussMap(r)
  r <- r/maxValue(r)
  abund <- raster(r)
  #abund[mid] <- 10000

  spiralsD <- spiralDistances(r, maxDis = 70000, cellSize = res(massAttacksMap)[1])
  tos <- ( spiralsD[, c("row", "col")]) * res(r)
  colnames(tos) <- c("x", "y")

  sdDist = p[4]
  par(mfrow = c(4,4))
  for (i in 1:4) {
    mid <- SpaDES.tools::middlePixel(abund)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    meanDist = rlnorm(1, log(p[[1]]), log(p[[5]]))
    asymDir <- rnorm(1, p["advectionDir"], p["advectionDirSD"]*3)
    dd <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                           meanDist = meanDist, sdDist = p["sdDist"], landscape = r,
                           propPineMapInner = r, fromCell = mid,
                           toCells = toCells,
                           asym = p["advectionMag"], angle = ddd[, "angles"],
                           asymDir = asymDir)
    abund[toCells] <- a1 <- -dd * 10000/(sum(-dd, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd1 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                            meanDist = meanDist, sdDist = p["sdDist"], landscape = r,
                            propPineMapInner = r, fromCell = mid,
                            toCells = toCells,
                            asym = p["advectionMag"], angle = ddd[, "angles"],
                            asymDir = asymDir)
    abund[toCells] <- a2 <- -dd1 * 10000/(sum(-dd1, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd2 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                            meanDist = meanDist, sdDist = p["sdDist"], landscape = r,
                            propPineMapInner = r, fromCell = mid,
                            toCells = toCells,
                            asym = p["advectionMag"], angle = ddd[, "angles"],
                            asymDir =  asymDir)
    abund[toCells] <-  a3 <- -dd2 * 10000/(sum(-dd2, na.rm = TRUE))
    plot(abund)
    a4 <- rescale(a1, c(0,1)) + rescale(a2, c(0,1)) + rescale(a3, c(0,1))
    abund[toCells] <- a4
    plot(abund)
  }
  plot(type = "n", x = 1:10, axes = F, xlab = "", ylab = "")
  legend(x = 2, y = 8, col = c("red", "green"), lty= 1,
         legend = c("annual average", "global average"),
         cex = 2)
  mtext(text = "Dispersal kernels", outer = TRUE, line = -3, xpd = TRUE)
  # HAVEN'T MADE IT PAST HERE

  ##########################################
  maxD <- 70000 + 125
  r <- raster(extent(-maxD+ 250, maxD + 250 , -maxD, maxD), res = 250)
  r[] <- 1
  r <- gaussMap(r)
  r <- r/maxValue(r)
  abund <- raster(r)
  #abund[mid] <- 10000

  spiralsD <- spiralDistances(r, maxDis = 70000, cellSize = res(massAttacksMap)[1])
  tos <- ( spiralsD[, c("row", "col")]) * res(r)
  colnames(tos) <- c("x", "y")

  sdDist = p[4]
  # p["advectionMag"] <- 15000
  par(mfrow = c(4,4))
  for (i in 1:4) {
    mid <- SpaDES.tools::middlePixel(abund)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    meanDist = rlnorm(1, log(p[[1]]), log(p[[5]]))
    asymDir <- rnorm(1, p["advectionDir"], p["advectionDirSD"]*3)
    dd <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                           meanDist = meanDist, sdDist = p["sdDist"], landscape = r,
                           propPineMapInner = r, fromCell = mid,
                           toCells = toCells,
                           asym = p["advectionMag"], angle = ddd[, "angles"],
                           asymDir = asymDir)
    abund[toCells] <- a1 <- -dd * 10000/(sum(-dd, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd1 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                           meanDist = meanDist, sdDist = p["sdDist"], landscape = r,
                           propPineMapInner = r, fromCell = mid,
                           toCells = toCells,
                           asym = p["advectionMag"], angle = ddd[, "angles"],
                           asymDir = asymDir)
    abund[toCells] <- a2 <- -dd1 * 10000/(sum(-dd1, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd2 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                           meanDist = meanDist, sdDist = p["sdDist"], landscape = r,
                           propPineMapInner = r, fromCell = mid,
                           toCells = toCells,
                           asym = p["advectionMag"], angle = ddd[, "angles"],
                           asymDir =  asymDir)
    abund[toCells] <-  a3 <- -dd2 * 10000/(sum(-dd2, na.rm = TRUE))
    plot(abund)
    a4 <- rescale(a1, c(0,1)) + rescale(a2, c(0,1)) + rescale(a3, c(0,1))
    abund[toCells] <- a4
    plot(abund)
  }


  ###############################################


  stop()

  # Visualize with maps
  # Kernel
  advectionMagTmp <- p[2] ;
  advectionDir <- p[3]
  r <- raster(extent(-20000, 20000, -20000, 20000), res = 100)
  r[] <- 1
  abund <- raster(r)
  mid <- SpaDES.tools::middlePixel(abund)
  abund[mid] <- 10000
  out <- spread3(mid, rasQuality = r, rasAbundance = abund,
          advectionDir = advectionDir,
          advectionMag = advectionMagTmp,
          meanDist = meanDist,
          sdDist = sdDist, dispersalKernel = "weibull", plot.it = F)
  abund[out$pixels] <- out$abundSettled
  abund[mid] <- max(abund[], na.rm = TRUE) * 1.5
  clearPlot(); Plot(abund)

  #########################################################
  #################### PREDICT MAPS
  #########################################################
  nams <- names(massAttacksMap)
  for (i in seq(nams)[-length(nams)]) {
    a <- massAttacksMap[[nams[i]]]
    b <- propPineMap
    whTos <- which(b[] >= 0.4)
    whFroms <- which(a[] > 0)
    if (length(whFroms) && length(whTos)) {
      env <- new.env()
      env$atksRasNextYr <- b
      env$atksKnownNextYr <- data.table(pixels = whTos)
      env$starts <- whFroms
      env$currentAttacks <- a
      env$dispersalKernel <- "exponential"
      env$propPineMapInner <- propPineMap
      env$p <- p
      st <- system.time(out <- eval(quotedSpread, envir = env))
      pixels <- cellFromXY(massAttacksMap, out[, c("x", "y")])
      predictedMap <- raster(massAttacksMap)
      predictedMap[whTos] <- log(-out[, "val"] + 1)
      clearPlot()
      currentAttacks <- env$currentAttacks
      nextYearAttacks <- massAttacksMap[[nams[i+1]]]
      Plot(currentAttacks, nextYearAttacks, predictedMap)
    }
  }
  #########################################################
  #########################################################


  nams <- names(massAttacksMap)
  # diffs <- list()
  diffs <- lapply(seq(nams)[-length(nams)], function(i) {
    a <- massAttacksMap[[nams[i]]]
    b <- propPineMap
    whTos <- which(b[] >= 0.2)
    whFroms <- which(a[] > 0)
    if (length(whFroms) && length(whTos)) {
      env <- new.env()
      env$atksRasNextYr <- b
      env$atksKnownNextYr <- data.table(pixels = whTos)
      env$starts <- whFroms
      env$currentAttacks <- a
      env$dispersalKernel <- "exponential"
      env$propPineMapInner <- propPineMap
      env$p <- p
      st <- system.time(out <- eval(quotedSpread, envir = env))
      pixels <- cellFromXY(massAttacksMap, out[, c("x", "y")])
      predictedMap <- raster(massAttacksMap)
      predictedMap[whTos] <- log(-out[, "val"] + 1)
      clearPlot()
      currentAttacks <- env$currentAttacks
      nextYearAttacks <- massAttacksMap[[nams[i+1]]]
      Plot(currentAttacks, nextYearAttacks, predictedMap)


      dirs <- distanceFromEachPoint(froms, tos, a, angles = TRUE)

      browser()
      pixels <- cellFromXY(a, xy = dirs[, c("x", "y")])
      dt <- data.table(pixels, dirs[, c("dists", "angles")])
      setorder(dt, dists)
      dirs <- dt[, .SD[1:10], by = "pixels"]

      # From https://www.themathdoctors.org/averaging-angles/
      avgDir <- try(CircStats::deg(atan(sum(sin(dirs[, "angles"]))/sum(cos(dirs[, "angles"])))))
      hist(dirs[["angles"]] %% (2*pi))
      avgDist <- mean(dirs$dists)
      if (is(avgDir, "try-error")) browser()
    } else {
      avgDir <- NA
      avgDist <- NA
    }

    a <- trim(a)
    a[] <- a[] > 0
    b <- trim(b)
    b[] <- (b[] > 0) + 1
    d <- mosaic(a, b, fun = sum)
    d[d[]==0] <- NA
    levels(d) <- data.frame(ID = 0:3, c("none", "prevYr", "nextYr", "both"))
    setColors(d) <- c("transparent", "red", "blue", "purple")
    d <- setMinMax(d)
    writeRaster(d, filename = paste0("MPB_Yrs_", nams[i], "-", nams[i+1], ".tif"), overwrite = TRUE)
    # png(filename = paste0("MPB_Yrs_", nams[i], "-", nams[i+1], ".png"),
    #     width = 2000)
    # Plot(d)
    # dev.off()
    # diffs[[i]] <- d
    list(ras = d, avgDir = avgDir, avgDist = avgDist)
  })

  diffs2 <- purrr::transpose(diffs)
  diffs <- diffs2$ras
  dirAvgs <- diffs2$avgDir
  distAvgs <- unlist(diffs2$avgDist)
  dirAvgs <- data.table(period = paste0(nams[-length(nams)], "-to-", nams[-1]),
                        dirAvgs, distAvgs)
  dirAvgs <- as.data.table(dirAvgs, keep.rownames = "YearPair")
  names(dirAvgs) <-
  names(diffs) <- nams[-length(nams)]
  dev()
  clearPlot()
  Plot(diffs, visualSqueeze = 1)

  # Average direction


  mn <- (meanDist + advectionMagTmp)
  sd <- mn/sdDist # 0.8 to 2.0 range
  shape <- (sd/mn)^(-1.086)
  scale <- mn/exp(lgamma(1+1/shape))
  aa <- curve(dweibull(x, shape = shape, scale = scale), from = 0, 1e5)

  mn <- max(10, (meanDist - advectionMagTmp))
  sd <- mn/sdDist # 0.8 to 2.0 range
  shape <- (sd/mn)^(-1.086)
  scale <- mn/exp(lgamma(1+1/shape))
  bb <- curve(dweibull(x, shape = shape, scale = scale), from = 0, 1e5)
  clearPlot()
  Plot(aa$x, aa$y, addTo = "Kernel_with_wind", type = "l")
  Plot(bb$x, bb$y, addTo = "Kernel_without_wind", type = "l")



  ## sum negative log likelihood for attacked pixels


  ## TODO: something other than simple sum of squares?
  #  metric <- (atkAreaData - atkAreaSim)^2 #+ (SNLL / 10^3)



  if (isTRUE(!is.na(params$.plotInitialTime))) {

    out2 <- out[, list(abundSettled = sum(abundSettled)), by = c("distance", "direction")]
    maxDist <- max(out$distance)
    rr <- raster::raster(extent(-maxDist, maxDist, -maxDist, maxDist), res = res(propPineMap))
    out3 <- out2[abundSettled > 0]
    xx = out3$distance * (sin(out3$direction))
    yy = out3$distance * (-cos(out3$direction))
    pixels <- cellFromXY(rr, cbind(xx, yy))
    rr[pixels] <- out3$abundSettled
    raster::plot(rr)


    out2 <- out[, list(abundSettled = sum(abundSettled)), by = c("pixels")]
    rr <- raster::raster(propPineMap)
    rr[out2$pixels] <- out2$abundSettled
    rr2 <- trim(rr)
    rr3 <- raster(rr2)
    rr3[rr2[] > 10] <- 2
    rr4 <- trim(rr3)
    rr5 <- crop(rr2, rr4)
    Plot(rr5)


  }
  if (EliotTesting) {
    fname <- file.path(outputPath(sim), paste0("spread", current(sim)$eventTime, ".tif"))
    tf <- raster::rasterTmpFile()
    r <- stack(saveStack)
    r2 <- calc(r, sum, filename = tf, overwrite = TRUE)
    writeRaster(r2, fname, overwrite = TRUE)
    unlink(saveStack, tf)
  }

  # if (EliotTesting) {
  #   tmpStackObj <- stack(saveStack)
  #   ex <- extent(tmpStackObj)
  #   ex@ymax <- ex@ymax - 4000
  #   ex@ymin <- 7389000
  #   ex@xmax <- -901000
  #   out2 <- crop(tmpStackObj, ex) * 10
  #   if (require(animation)) {
  #     gifName <- "C:\\Eliot\\Google Drive\\McIntire-lab\\figures\\MPB animation 700px small area 3.gif"
  #     #gifName <- file.path(tempdir(), "animation.gif")
  #     saveGIF(interval = 0.1, ani.height = 700, ani.width = 700, movie.name = gifName, expr = {
  #       for (i in seq(numLayers(out2))) plot(out2[[i]])
  #     })
  #   }
  #   stop("End it here")
  # }

  migrantsDT <- out[, list(NEWATKs = sum(abundSettled)), by = "pixels"]
  massAttacksDTNew <- massAttacksDT[migrantsDT, on = c(ID = "pixels")]

  # ! ----- STOP EDITING ----- ! #
  return(invisible(massAttacksDTNew))
}



objFun <- function(p, propPineMap, currentAttacks, advectionDir, advectionMag,
                   minNumAgents, massAttacksMap, currentTime, reps = 10,
                   quotedSpread, fitType = "ss1", omitPastPines, sdDist, dispersalKernel) {

  # DEoptim, when used across a network, inefficiently moves large objects every time
  #   it calls makes an objFun call. Better to move the objects to each node's disk,
  #   then read them in locally from disk --> faster when bandwidth is slow
  if (missing(propPineMap))
    propPineMap <- get("propPineMap", envir = .GlobalEnv)
  if (missing(massAttacksMap))
    massAttacksMap <- get("massAttacksMap", envir = .GlobalEnv)
  if (missing(currentAttacks))
    currentAttacks <- get("currentAttacks", envir = .GlobalEnv)
  if (missing(advectionDir))
    advectionDir <- get("advectionDir", envir = .GlobalEnv)
  if (missing(advectionMag))
    advectionMag <- get("advectionMag", envir = .GlobalEnv)
  if (missing(minNumAgents))
    minNumAgents <- get("minNumAgents", envir = .GlobalEnv)
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

  if (fitType == "likelihood") {
    if (reps < 10) {
      message("reps must be at least 10 if using likelihood; setting to 10")
      reps <- 10
    }
  }

  allYears <- names(massAttacksMap)

  combinations <- seq_along(allYears[-1])
  years <- data.frame(
    startYears = allYears[-length(allYears)],
    endYears = allYears[-1])
  # combinations <- cbind(combinations, years[combinations[, "years"], ])
  # combinations <- as.data.frame(combinations[, c("reps", "startYears", "endYears")])

  maxObjFunValIndiv <- -Inf
  # mam <- raster::unstack(massAttacksMap)
  outBig <- purrr::pmap(
    .l = years,
    starts = starts,
              #currentAttacks = mam[-length(mam)],
              #atksRasNextYr = mam[-1]
    #startYears = startYears,
    #endYears = endYears,
    #starts = startsOuter,
    massAttacksMap = massAttacksMap,
    #advDir = advectionDir,
    #advMag = advectionMag,
    p = p, reps = reps,
    minNumAgents = minNumAgents,
    propPineMapInner = propPineMap,
    quotedSpread = quotedSpread,
    fitType = fitType,
    omitPastPines = omitPastPines,
    .f = objFunInner
  )

  sum(unlist(outBig), na.rm = TRUE)

}

objFunInner <- function(reps, starts, startYears, endYears, p, minNumAgents,
                        # atksRasNextYr, currentAttacks,
                        massAttacksMap,
                        propPineMapInner, objFunValOnFail = 1e3,
                        # advDir, advMag,
                        quotedSpread, fitType, omitPastPines) {
  currentAttacks <- massAttacksMap[[startYears]]

  starts <- which(massAttacksMap[[startYears]][] > 0)

  if (omitPastPines) {

    namesMAM <- names(massAttacksMap)

    namesToStartYears <- namesMAM[seq(which(namesMAM %in% startYears))]
    # unlist them
    mam <- if (length(namesToStartYears) == 1)
      list(massAttacksMap[[namesToStartYears]])
    else
      raster::unstack(massAttacksMap[[namesToStartYears]])
    names(mam) <- namesToStartYears
    massAttacksDTYearsToPres <- lapply(mam, function(x) {
      pixels = which(x[] > 0)
      setDT(list(pixels = pixels, ATKTREES = x[][pixels]))
    })
    massAttacksDTYearsToPres <- rbindlist(massAttacksDTYearsToPres, idcol = "Year")
  }
  atksRasNextYr <- massAttacksMap[[endYears]]
  wh <- which(atksRasNextYr[] > 0)
  atksKnownNextYr <- setDT(list(pixels = wh, ATKTREES = atksRasNextYr[][wh]))

  env <- environment()

  out <- lapply(seq_len(reps), function(rep) as.data.table(eval(quotedSpread, envir = env)))
  out <- rbindlist(out, idcol = "rep")

  if (fitType == "distanceFromEachPoint") {
    objFunVal <- out[,"val"]

  } else {



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
  objFunVal <- round(sum(objFunVal, na.rm = TRUE), 3)
  print(paste("Done startYear: ", startYears, " ObjFunVal: ", objFunVal))
  return(objFunVal)
}

spiralDistances <- function(pixelGroupMap, maxDis, cellSize) {
  spiral <- which(focalWeight(pixelGroupMap, maxDis, type = "circle")>0, arr.ind = TRUE) -
    ceiling(maxDis/cellSize) - 1
  spiral <- cbind(spiral, dists = sqrt( (0 - spiral[,1]) ^ 2 + (0 - spiral[, 2]) ^ 2))
  spiral <- spiral[order(spiral[, "dists"], apply(abs(spiral), 1, sum),
                         abs(spiral[, 1]), abs(spiral[, 2])),, drop = FALSE]
}

distanceFunction <- function(dist, angle, landscape, fromCell, toCells, nextYrVec,
                             propPineMapInner, asym, asymDir, meanDist, sdDist, dispersalKernel) {

  if (asym > 0) {
    fromXY <- xyFromCell(propPineMapInner, fromCell)
    toXYs <-  xyFromCell(propPineMapInner, toCells)
    advectionXY <- c(x = sin(rad(asymDir))*asym, y = cos(rad(asymDir))*asym)
    # neededXY <- toXYs - rep(advectionXY, each = NROW(toXYs))

    # x <- numeric(NROW(toXYs))
    # y <- x
    x <- toXYs[, 'x'] - advectionXY[1]/meanDist*dist
    y <- toXYs[, 'y'] - advectionXY[2]/meanDist*dist
    newTo <- cbind(x, y)
    colnames(newTo) <- c("x", "y")
    # xLT0 <- toXYs[, "x"] < 0
    # whXLT0 <- which(xLT0)
    # whXGT0 <- which(!xLT0)
    # x[whXLT0] <- toXYs[whXLT0, 'x'] - advectionXY[1]/meanDist*dist[whXLT0]
    # x[whXGT0] <- toXYs[whXGT0, 'x'] - advectionXY[1]/meanDist*dist[whXGT0]
    # yLT0 <- toXYs[, "y"] < 0
    # whYLT0 <- which(yLT0)
    # whYGT0 <- which(!yLT0)
    # y[whYLT0] <- toXYs[whYLT0, 'y'] - advectionXY[2]/meanDist*dist[whYLT0]
    # y[whYGT0] <- toXYs[whYGT0, 'y'] - advectionXY[2]/meanDist*dist[whYGT0]
    # x <- toXYs[, 'x'] - advectionXY[1]/meanDist*toXYs[, 'x']
    # y <- toXYs[, 'y'] - advectionXY[2]/meanDist*toXYs[, 'y']
    dists <- cbind(origX = toXYs[, "x"], origY = toXYs[, "y"], origDist = dist,
                   .pointDistance(from = fromXY, newTo, angles = FALSE))
    dist <- dists[, "dists"]
  }

  if (grepl("Weibull", dispersalKernel)) {
    mn <- (meanDist)
    sd <- mn/sdDist # 0.8 to 2.0 range
    shape <- (sd/mn)^(-1.086)
    scale <- mn/exp(lgamma(1+1/shape))

    if (grepl("Weibull", dispersalKernel)) {
      prob <- -dweibull(dist, scale = scale, shape = shape)
    }

    infs <- is.infinite(prob)
    if (any(infs))
      prob[infs] <- 0
  } else {
    prob <- -dexp(dist, rate = 1/meanDist)

  }

  prob <- prob * landscape[fromCell] * propPineMapInner[][toCells]
  # if (any(is.infinite(out1))) browser()
  # angle <- -30:30/10
  # plot(deg(angle), asym^sin(angle-rad(asymDir)))
  # plot(deg(angle), asym^cos(angle-rad(asymDir)))
}
