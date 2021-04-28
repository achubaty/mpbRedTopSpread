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
  reqdPkgs = list("achubaty/amc@development", "data.table", "quickPlot",
                  "PredictiveEcology/LandR@development",
                  "raster", "RColorBrewer", "reproducible"),
  parameters = rbind(
    defineParameter("advectionDir", "numeric", 90, 0, 359.9999,
                    "The direction of the spread bias, in degrees from north"),
    defineParameter("advectionMag", "numeric", 1000, NA, NA,
                    "The magnitude of the directional bias of spread"),
    defineParameter("bgSettlingProp", "numeric", 0.1, 0, 1,
                    "The proportion of beetles that settle from those that could potentially settle, even if no pine"),
    defineParameter("meanDist", "numeric", 1000, NA, NA,
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
                                     currentAttacks = sim$currentAttacks,
                                     params = P(sim)
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
dispersal <- function(sim) {
  ## check that MPB and pine rasters are the same resolution and ncells
  if (fromDisk(sim$pineMap))
    sim$pineMap[] <- sim$pineMap[]

  ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (2017)
  MAXTREES <- 1125 * prod(res(sim$pineMap)) / 100^2 ## TODO: round this?

  ## asymmetric spread (biased eastward)
  # lodgepole pine and jack pine together
  propPineMap <- if (nlayers(sim$pineMap) > 1) sum(sim$pineMap) else sim$pineMap
  # mv <- maxValue(propPineMap)
  #
  # if (mv > 1)
  #   propPineMap[] <- propPineMap[] / 100 ## much faster than calc; ## TODO: allow different params per species

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
  #propPineMap[] <- pmin(1, propPineMap[] + P(sim)$bgSettlingProp) ## TODO: why divide by 100??

  if (EliotTesting) { # TODO -- delete EliotTesting when no longer desired
    a <- extent(sim$studyArea)
    starts <- sim$massAttacksDT[ATKTREES > 0]$ID
    d <- raster(sim$currentAttacks)
    d[starts] <- 1
    d <- crop(d, a)
    starts <- which(d[] > 0)
    currentAttacks <- crop(sim$currentAttacks, a) * 300
    propPineMap <- crop(propPineMap, a)
    saveStack <- raster::rasterTmpFile()
  } else {
    minNumAgents <- 50
    starts <- sim$massAttacksDT[["ID"]][sim$massAttacksDT$ATKTREES > 0]
    saveStack <- NULL
    currentAttacks <- sim$currentAttacks
  }

  #  st1 <- system.time({
  out <- SpaDES.tools::spread3(start = starts,
                               rasQuality = propPineMap,
                               rasAbundance = currentAttacks,
                               advectionDir = P(sim)$advectionDir,
                               advectionMag = P(sim)$advectionMag,
                               meanDist = P(sim)$meanDist,
                               plot.it = !is.na(P(sim)$.plotInitialTime),
                               minNumAgents = minNumAgents,
                               verbose = 2,
                               skipChecks = TRUE,
                               saveStack = saveStack) ## saveStack is the filename to save to
  #})
  if (EliotTesting) {
    fname <- file.path(outputPath(sim), paste0("spread", current(sim)$eventTime, ".tif"))
    tf <- raster::rasterTmpFile()
    r <- stack(saveStack)
    r2 <- calc(r, sum, filename = tf, overwrite = TRUE)
    writeRaster(r2, fname, overwrite = TRUE)
    unlink(saveStack, tf)
  }

  if (FALSE) {
    atks <- simOut$massAttacksDT
    nPix <- atks[ATKTREES > 0, .N]
    atkAreaSim <- nPix * prod(res(simOut$rasterToMatch)) / (100^2) ## area in ha

    ## attacked area from data
    atksRas <- simOut$massAttacksMap[[paste0("X", timesFit$end)]]
    atks <- data.table(ID = 1L:ncell(atksRas), ATKTREES = atksRas[])
    nPix <- atks[ATKTREES > 0, .N] ## total number of pixels
    atkAreaData <- nPix * prod(res(simOut$rasterToMatch)) / (100^2) ## area in ha

    ## sum negative log likelihood for attacked pixels


    ## TODO: something other than simple sum of squares?
    metric <- (atkAreaData - atkAreaSim)^2 #+ (SNLL / 10^3)
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
  out2 <- sim$massAttacksDT[migrantsDT, on = c(ID = "pixels")]

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

dispersal2 <- function(pineMap, studyArea, massAttacksDT, currentAttacks, params) {
  ## check that MPB and pine rasters are the same resolution and ncells
  if (fromDisk(pineMap))
    pineMap[] <- pineMap[]

  ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (2017)
  MAXTREES <- 1125 * prod(res(pineMap)) / 100^2 ## TODO: round this?

  ## asymmetric spread (biased eastward)
  # lodgepole pine and jack pine together
  propPineMap <- if (nlayers(pineMap) > 1) sum(pineMap) else pineMap
  # mv <- maxValue(propPineMap)
  #
  # if (mv > 1)
  #   propPineMap[] <- propPineMap[] / 100 ## much faster than calc; ## TODO: allow different params per species

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
  #propPineMap[] <- pmin(1, propPineMap[] + P(sim)$bgSettlingProp) ## TODO: why divide by 100??

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

  #  st1 <- system.time({
  out <- SpaDES.tools::spread3(start = starts,
                               rasQuality = propPineMap,
                               rasAbundance = currentAttacks,
                               advectionDir = params$advectionDir,
                               advectionMag = params$advectionMag,
                               meanDist = params$meanDist,
                               plot.it = !is.na(params$.plotInitialTime),
                               minNumAgents = minNumAgents,
                               verbose = 2,
                               skipChecks = TRUE,
                               saveStack = saveStack) ## saveStack is the filename to save to
  #})
  if (EliotTesting) {
    fname <- file.path(outputPath(sim), paste0("spread", current(sim)$eventTime, ".tif"))
    tf <- raster::rasterTmpFile()
    r <- stack(saveStack)
    r2 <- calc(r, sum, filename = tf, overwrite = TRUE)
    writeRaster(r2, fname, overwrite = TRUE)
    unlink(saveStack, tf)
  }

  if (FALSE) {
    atks <- simOut$massAttacksDT
    nPix <- atks[ATKTREES > 0, .N]
    atkAreaSim <- nPix * prod(res(simOut$rasterToMatch)) / (100^2) ## area in ha

    ## attacked area from data
    atksRas <- simOut$massAttacksMap[[paste0("X", timesFit$end)]]
    atks <- data.table(ID = 1L:ncell(atksRas), ATKTREES = atksRas[])
    nPix <- atks[ATKTREES > 0, .N] ## total number of pixels
    atkAreaData <- nPix * prod(res(simOut$rasterToMatch)) / (100^2) ## area in ha

    ## sum negative log likelihood for attacked pixels


    ## TODO: something other than simple sum of squares?
    metric <- (atkAreaData - atkAreaSim)^2 #+ (SNLL / 10^3)
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
  out2 <- massAttacksDT[migrantsDT, on = c(ID = "pixels")]

  # ! ----- STOP EDITING ----- ! #
  return(invisible())
}
