defineModule(sim, list(
  name = "mpbRedTopSpread",
  description = "Mountain Pine Beetle Red Top Growth Model: Short-run Potential for Establishment, Eruption, and Spread",
  keywords = c("mountain pine beetle, outbreak dynamics, eruptive potential, spread, climate change, twitch response"),
  authors = c(
    person(c("Alex", "M"), "Chubaty", email = "alexander.chubaty@canada.ca", role = c("aut", "cre")),
    person(c("Eliot", "J B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut")),
    person(c("Barry", "J"), "Cooke", email = "barry.cooke@canada.ca", role = c("aut"))
  ),
  childModules = character(),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list(),
  reqdPkgs = list("amc", "data.table", "quickPlot", "raster", "RColorBrewer", "reproducible"),
  parameters = rbind(
    defineParameter("advectionDir", "numeric", 90, 0, 359.9999, "The direction of the spread bias, in degrees from north"),
    defineParameter("advectionMag", "numeric", 3000, NA, NA, "The magnitude of the directional bias of spread"),
    defineParameter("bgSettlingProp", "numeric", 0.1, 0, 1, "The proportion of beetles that settle from those that could potentially settle, even if no pine"),
    defineParameter("meanDist", "numeric", 1000, NA, NA, "Expected dispersal distance (m); ~63% go less than this distance"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated?")
  ),
  inputObjects = bind_rows(
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
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureStandVolume.tar")
  ),
  outputObjects = bind_rows(
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
      sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "dispersal",
                           eventPriority = 4.5)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbRedTopSpread", "save")
    },
    "dispersal" = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- dispersal(sim)

      sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "dispersal",
                           eventPriority = 4.5)

      # ! ----- STOP EDITING ----- ! #
    },
    "plot" = {
      # ! ----- EDIT BELOW ----- ! #

      sim <- plotFn(sim)

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  return(invisible(sim))
}

### initilization
Init <- function(sim) {
  ## stand age map
  if (!suppliedElsewhere("standAgeMap", sim)) {
    standAgeMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.tif")
    sim$standAgeMap <- Cache(prepInputs,
                             targetFile = basename(standAgeMapFilename),
                             archive = asPath(c("kNN-StructureStandVolume.tar",
                                                "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.zip")),
                             destinationPath = dPath,
                             url = na.omit(extractURL("standAgeMap")),
                             fun = "raster::raster",
                             studyArea = sim$studyArea,
                             #rasterToMatch = sim$rasterToMatch,
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = paste0(tools::file_path_sans_ext(basename(standAgeMapFilename)), "_cropped"),
                             overwrite = TRUE,
                             userTags = c("stable", currentModule(sim)))
    sim$standAgeMap[] <- asInteger(sim$standAgeMap[])
  }

  ## raster to match
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- sim$standAgeMap
  }

  return(invisible(sim))
}

### plotting
plotFn <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event

  #Plot(amc::dt2raster(sim$massAttacksDT, sim$massAttacksMap, "ATKTREES")) ## TODO: fix this

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### spread
dispersal <- function(sim) {
  ## check that MPB and pine rasters are the same resolution and ncells
  stopifnot(all.equal(res(sim$currentAttacks), res(sim$pineMap))) # TODO: use rasterToMatch

  ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (unpublished)
  MAXTREES <- 1125 * prod(res(sim$pineMap)) / 100^2 ## TODO: round this?

  ## asymmetric spread (biased eastward)
  # lodgepole pine and jack pine together ## TODO: allow different parameterizations per species
  propPineMap <- sim$pineMap[["Pinu_sp"]]# / 100
  propPineMap[is.na(propPineMap[])] <- 0
  if (exists("EliotTesting")) {
    EliotTesting <- TRUE
    sim@params$mpbRedTopSpread$bgSettlingProp <- 0.4
    sim@params$mpbRedTopSpread$advectionMag <- 30000
    sim@params$mpbRedTopSpread$.plotInitialTime <- NA
  } else {
    EliotTesting <- FALSE
  }
  propPineMap[] <- pmin(1, propPineMap[]/100 + P(sim)$bgSettlingProp)

  if (EliotTesting) { # TODO -- delete EliotTesting when no longer desired
    a <- extent(sim$studyArea)
    starts <- sim$massAttacksDT[ATKTREES > 0]$ID
    d <- raster(sim$currentAttacks)
    d[starts] <- 1
    d <- crop(d, a)
    starts <- which(d[] > 0)
    b <- crop(sim$currentAttacks, a) * 100
    propPineMap <- crop(propPineMap, a)
    saveStack <- raster::rasterTmpFile()

  } else {
    starts <- sim$massAttacksDT[ATKTREES > 0]$ID
    saveStack <- NULL
    b <- sim$currentAttacks
  }

  st1 <- system.time(out <- SpaDES.tools::spread3(start = starts, ## TODO: remove the small subset 1:100
                        rasQuality = propPineMap,
                        rasAbundance = b,#sim$currentAttacks,
                        advectionDir = P(sim)$advectionDir,
                        advectionMag = P(sim)$advectionMag,
                        meanDist = P(sim)$meanDist,
                        plot.it = !is.na(P(sim)$.plotInitialTime),
                        minNumAgents = 100,
                        verbose = 2,
                        saveStack = saveStack)) ## saveStack is the filename to save to
  if (EliotTesting) {
    tmpStackObj <- stack(saveStack)
    ex <- extent(tmpStackObj)
    ex@ymax <- ex@ymax - 1000
    ex@ymin <- 7389000
    ex@xmax <- -901000
    out2 <- crop(tmpStackObj, ex)
    if (require(animation)) {
      gifName <- file.path(tempdir(), "animation.gif")
      saveGIF(interval = 0.1, ani.height = 700, ani.width = 700, movie.name = gifName, expr = {
        for (i in seq(numLayers(out2))) plot(out2[[i]])
      })
    }
    stop("End it here")
  }
  migrantsDT <- out[, list(NEWATKs = sum(abundSettled)), by = "pixels"]
  out2 <- sim$massAttacksDT[migrantsDT, on = c(ID = "pixels")]


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

