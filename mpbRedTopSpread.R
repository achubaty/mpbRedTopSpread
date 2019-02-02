
# Everything in this file gets sourced during simInit, and all functions and objects
#  are put into the simList. To use objects and functions, use sim$xxx.
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
    defineParameter("asymmetry", "numeric", 2, NA, NA, "The magnitude of the directional bias of spread"),
    defineParameter("asymmetryAngle", "numeric", 90, NA, NA, "The direction of the spread bias, in degrees from north"),
    defineParameter("dispersalInterval", "numeric", 1, NA, NA, "This describes the interval time between dispersal events"),
    defineParameter("dispersalKernel", "character", "NegExp", NA, NA, "Name of the dispersal kernel to use"),
    defineParameter("dispersalKernelLambda", "numeric", 1, NA, NA, "Dispersal kernel lambda parameter"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated?")
  ),
  inputObjects = bind_rows(
    expectsInput("massAttacksDT", "data.table",
                 desc = "Current MPB attack map (number of red attacked trees).",
                 sourceURL = NA),
    expectsInput("massAttacksMap", "RasterStack",
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
      sim <- scheduleEvent(sim, time(sim) + P(sim)$dispersalInterval,
                           "mpbRedTopSpread", "dispersal")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbRedTopSpread", "save")
    },
    "dispersal" = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      sim <- dispersal(sim)

      sim <- scheduleEvent(sim, time(sim) + P(sim)$dispersalInterval,
                           "mpbRedTopSpread", "dispersal")

      # ! ----- STOP EDITING ----- ! #
    },
    "plot" = {
      # ! ----- EDIT BELOW ----- ! #

      plot(amc::dt2raster(sim$mpbSpreadDT))

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
  ## dispersal kernels
  mod$dispKern <- switch(
    P(sim)$dispersalKernel,
    "NegExp" = function(disFar, disNear, lambda) {
      (1 - exp(-lambda * disFar)) - (1 - exp(-lambda * disNear))
    }
  )

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
  #Plot("object")

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### spread
dispersal <- function(sim) {
  ## check that MPB and pine rasters are the same resolution and ncells

  stopifnot(all.equal(res(sim$massAttacksMap), res(sim$pineMap))) # TODO: use rasterToMatch

  ## X and Y map resolutions should be equal (square pixels)
  MAPRES <- unique(res(sim$massAttacksMap))
  stopifnot(length(MAPRES) == 1)

  ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (unpublished)
  MAXTREES <- round(1125 * (MAPRES / 100) ^ 2)

  ASYMM <- 5     ## TODO: parameterize this
  BIAS <- 90     ## TODO: parameterize this
  LAMBDA <- 0.12 ## TODO: parameterize this
  SATDENS <- 332 ## TODO: parameterize this

  DEBUG <- FALSE ## TODO: parameterize this

  ## asymmetric spread (biased eastward)
  # lodgepole pine and jack pine together ## TODO: allow different parameterizations per species
  propPineMap <- (sim$pineMap[["Jack_Pine"]] + sim$pineMap[["Lodgepole_Pine"]]) / 100
  insect_spread(r = propPineMap, loci = sim$massAttacksDT$ID,
                asymmetry = ASYMM, asymmetryAngle = BIAS, lambda = LAMBDA,
                saturationDensity = SATDENS, total = MAXTREES, Ncpus = 2, debug = DEBUG)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
