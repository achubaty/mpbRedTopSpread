## package BioSIM not available on CRAN nor GitHub; needs to be installed as follows:
## install.packages("https://sourceforge.net/projects/repiceasource/files/latest", repos = NULL,  type = "source")
## install.packages("https://sourceforge.net/projects/biosimclient.mrnfforesttools.p/files/latest", repos = NULL,  type = "source")

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
  reqdPkgs = list("achubaty/amc@development",
                  "BioSIM", ##  ## not on CRAN/GitHub; see install info at top of this file
                  "CircStats", "data.table", "DEoptim", "EnvStats",
                  "SpaDES.core (>= 1.0.8)",
                  "PredictiveEcology/LandR@LCC2010 (>= 1.0.4)",
                  "parallelly",
                  "PredictiveEcology/pemisc@development (>= 0.0.3.9001)",
                  "PredictiveEcology/mpbutils (>= 0.1.2)", "purrr",
                  "quickPlot", "raster", "RColorBrewer", "PredictiveEcology/reproducible",
                  "PredictiveEcology/SpaDES.tools@development (>= 0.3.7.9021)"),
  parameters = rbind(
    defineParameter("dataset", "character", "Boone2011", NA, NA, "Which dataset to use for stand dynamic model fitting. One of 'Boone2011' (default), 'Berryman1979_fit', or 'Berryman1979_forced'. Others to be implemented later."),
    defineParameter("growthInterval", "numeric", 1, NA, NA, "This describes the interval time between growth events"),
    defineParameter("advectionDir", "numeric", 90, 0, 359.9999,
                    "The direction of the spread bias, in degrees from north"),
    defineParameter("advectionMag", "numeric", 1000, NA, NA,
                    "The magnitude of the directional bias of spread"),
    defineParameter("bgSettlingProp", "numeric", 0.1, 0, 1,
                    "The proportion of beetles that settle from those that could potentially settle, even if no pine"),
    defineParameter("meanDist", "numeric", 1000, NA, NA,
                    "Expected dispersal distance (m); ~63% go less than this distance"),
    defineParameter("quotedSpread", "language", NULL, NA, NA,
                    "A spread function that is contained within a quote(). If left at NULL or NA, it will use the module version"),
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
    expectsInput("climateSuitabilityMaps", "RasterStack",
                 "A time series of climatic suitablity RasterLayers, each with previx 'X' and the year, e.g., 'X2010'"),
    # expectsInput("massAttacksRas", "RasterLayer",
    #              desc = "Current MPB attack map (number of red attacked trees).",
    #              sourceURL = NA),
    expectsInput("massAttacksStack", "RasterStack",
                 desc = "Historical MPB attack maps (number of red attacked trees).",
                 sourceURL = NA),
    expectsInput("massAttacksDT", "data.table",
                 desc = "Current MPB attack map (number of red attacked trees).",
                 sourceURL = NA),
    expectsInput("pineMap", "RasterLayer", "Percent cover maps by species (lodgepole and jack pine)."),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = "if not supplied, will default to standAgeMap", # TODO: description needed
                 sourceURL = NA),
    expectsInput("standAgeMap", "RasterLayer",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("windDirStack", "RasterStack",
                 desc = "RasterStack of wind maps for every location in the study area",
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("massAttacksDT", "data.table", "Current MPB attack map (number of red attacked trees)."),
    createsOutput("propPineRas", "RasterLayer",
                  paste("Proportion (not percent) cover map of *all* pine. This will",
                  " be the pixel-based sum if there are more than one layer")),
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
           sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "grow", eventPriority = 4.5)
           sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "dispersal", eventPriority = 4.5)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "plot", .last() - 1)
           sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbRedTopSpread", "save", .last())
         },
         "grow" = {
           # do stuff for this event
           sim <- grow(sim)

           # schedule future event(s)
           ## growth needs to happen after spread:
           sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "grow", eventPriority = 5.5)
         },
         "dispersal" = {
           out <- dispersal2(quotedSpread <- P(sim)$quotedSpread,
                             pineMap = sim$pineMap, studyArea = sim$studyArea,
                             massAttacksDT = sim$massAttacksDT,
                             massAttacksStack = sim$massAttacksStack,
                             # massAttacksRas = sim$massAttacksRas,
                             params = P(sim),
                             windDirStack = sim$windDirStack,
                             windSpeedStack = sim$windSpeedStack,
                             rasterToMatch = sim$rasterToMatch,
                             bgSettlingProp = P(sim)$bgSettlingProp,
                             type = P(sim)$type,
                             currentTime = time(sim),
                             reqdPkgs = reqdPkgs(module = currentModule(sim), modulePath = modulePath(sim))[[currentModule(sim)]]
           )
           sim$plotStack <- out
           # sim <- dispersal(sim)
           sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "dispersal", eventPriority = 5.5)
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
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  cPath <- cachePath(sim)
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  #mod$prj <- paste("+proj=aea +lat_1=47.5 +lat_2=54.5 +lat_0=0 +lon_0=-113",
  #                 "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  mod$prj <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                   "+x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

  ## load study area
  if (!suppliedElsewhere("studyArea")) {
    sim$studyArea <- mpbStudyArea(ecoregions = c(112, 120, 122, 124, 126), mod$prj,
                                  cPath, dPath)
  }

  ## raster to match
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- Cache(
      LandR::prepInputsLCC,
      year = 2005,
      destinationPath = dPath,
      studyArea = sf::as_Spatial(sim$studyArea)
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

  if (!suppliedElsewhere("windDirStack", sim)) {
    aggFact <- P(sim)$aggFact

    ## make coarser
    aggRTM <- raster::raster(sim$rasterToMatch)
    aggRTM <- raster::aggregate(aggRTM, fact = aggFact)
    aggRTM <- LandR::aggregateRasByDT(sim$rasterToMatch, aggRTM, fn = mean)

    dem <- Cache(prepInputsCanDEM,
                 studyArea = sim$studyArea,
                 rasterToMatch = aggRTM,
                 destinationPath = inputPath(sim))

    ## TODO: using "ClimaticWind_Annual"; use mean wind for summer [flight] months only
    windStk <- LandR::BioSIM_getWind(
      dem = dem,
      years = P(sim)$years,
      climModel = P(sim)$climateModel,
      rcp = P(sim)$climateScenario
    )

    windDirStack <- disaggregate(windStk, fact = aggFact)
    sim$windDirStack <- raster::stack(crop(windDirStack, sim$rasterToMatch))  ## TODO: speed and dir??

    if (!compareRaster(sim$windDirStack, sim$rasterToMatch, stopiffalse = FALSE)) {
      warning("wind raster is not same resolution as sim$rasterToMatch; please debug")
      browser() ## TODO: remove
    }

  }

  return(invisible(sim))
}

### initilization
Init <- function(sim) {
  cr <- compareRaster(sim$rasterToMatch, sim$massAttacksStack, sim$pineMap, orig = TRUE)
  if (!isTRUE(cr))
    stop("Rasters sim$rasterToMatch, sim$massAttacksStack, sim$pineMap must all have same metadata; they do not")


  # AGGREGATE FROM HERE
  if (fromDisk(sim$pineMap))
    sim$pineMap[] <- sim$pineMap[]

  ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (2017)
  MAXTREES <- 1125 * prod(res(sim$pineMap)) / 100^2 ## TODO: round this?

  ## asymmetric spread (biased eastward)
  # lodgepole pine and jack pine together
  propPineRas <- if (nlayers(sim$pineMap) > 1) sum(sim$pineMap) else sim$pineMap
  mv <- maxValue(propPineRas)

  if (mv > 1)
    propPineRas[] <- propPineRas[] / 100 ## much faster than calc; ## TODO: allow different params per species

  nas <- is.na(propPineRas[]) & !is.na(sim$rasterToMatch[])
  propPineRas[nas] <- 0


  nams <- names(sim$massAttacksStack)
  rasTmplate <- raster(sim$massAttacksStack)
  rasCoarse <- raster(raster::aggregate(rasTmplate, fact = 4))
  massAttacksList <- raster::unstack(sim$massAttacksStack)
  opts <- options(reproducible.cacheSpeed = "fast")
  on.exit(options(opts))
  mams <- raster::stack(
    Cache(lapply, massAttacksList, function(mam)
      aggregateRasByDT(mam, rasCoarse, fn = sum))
  )
  names(mams) <- names(sim$massAttacksStack)
  sim$massAttacksStack <- mams
  sim$propPineRas <- Cache(aggregateRasByDT, propPineRas, rasCoarse, fn = mean)
  sim$windDirStack <- Cache(aggregate, sim$windDirStack, res(rasCoarse)[1]/res(sim$windDirStack)[1])
  sim$windSpeedStack <- Cache(aggregate, sim$windSpeedStack, res(rasCoarse)[1]/res(sim$windSpeedStack)[1])
  resRatio <- res(rasCoarse)[1]/res(sim$climateSuitabilityMaps)[1]
  if (resRatio < 1) {
    climateSuitabilityMaps <- Cache(disaggregate, sim$climateSuitabilityMaps, fact = 10)
    climateSuitabilityMaps <- raster::stack(crop(climateSuitabilityMaps, sim$massAttacksStack))
    sim$climateSuitabilityMaps <- climateSuitabilityMaps
  }
  if (!compareRaster(sim$massAttacksStack,
                sim$propPineRas,
                sim$windDirStack,
                sim$windSpeedStack,
                sim$climateSuitabilityMaps))
    stop("The coarser resolution files need to be all the same resolution")

  options(opts)
  ## create a data.table consisting of the reduced map of current MPB distribution,
  ## proportion pine, and climatic suitability;
  ## use only the start year's non-zero and non-NA data
  ids <- sim$massAttacksStack[]
  ids <- which(ids > 0, arr.ind = TRUE)
  mpb.sp <- raster::xyFromCell(sim$massAttacksStack, cell = ids[, "row"])
  # crs(mpb.sp) <- crs(sim$massAttacksStack)
  sim$massAttacksDT <- as.data.table(ids)
  setnames(sim$massAttacksDT, old = c("row", "col"), new = c("pixel", "layerNum"))
  sim$massAttacksDT <- data.table(sim$massAttacksDT, mpb.sp)

  sim$massAttacksDT[, layerName := names(sim$massAttacksStack)[layerNum]]
  sim$massAttacksDT[, `:=`(
    CLIMATE = raster::extract(sim$climateSuitabilityMaps[[.BY[[1]]]], cbind(x, y)),
    ATKTREES = raster::extract(sim$massAttacksStack[[.BY[[1]]]], cbind(x, y))#,
    #windSpeed = raster::extract(sim$windSpeedStack[[.BY[[1]]]], cbind(x, y)),
    #windDir = raster::extract(sim$windDirStack[[.BY[[1]]]], cbind(x, y)),
    #propPineRas = raster::extract(sim$propPineRas, cbind(x, y))
    ),
    by = "layerName"]

  # ids <- which(!is.na(sim$currentAttacks[]) | (sim$currentAttacks[] > 0))
  # mpb.sp <- raster::xyFromCell(sim$currentAttacks, cell = ids)
# sim$massAttacksDT <- data.table(
  #   ID = ids[, "row"],
  #   layer = ids[, "col"],
  #   #X = mpb.sp[, 1],
  #   #Y = mpb.sp[, 2],
  #   ATKTREES = sim$massAttacksStack[ids[, "row"]],
  #   CLIMATE = raster::extract(sim$climateSuitabilityMaps, mpb.sp)
  # )

  ## growth data
  mod$growthData <- switch(P(sim)$dataset,
                           "Berryman1979_fit" = {
                             ## Berryman1979_forced
                             data.frame(
                               year = c(1:15),
                               log10Xtm1 = c(-3.1, -2.75, -2.7, -2.4, -2.3, -1.2, -1, 0.2, 0.9, 0.65,
                                             1.05, 0.95, 1.1, 1.5, 1.85),
                               log10Rt = c(0.35, 0.4, 0.1, -0.4, -0.65, 0.3, 1, 0.75, 1.2, -0.7,
                                           -0.4, 0.2, 0.45, 0.3, -0.78),
                               study = c(rep("Tunnock 1970", 9), rep("Parker 1973", 6)),
                               stringsAsFactors = TRUE
                             )
                           },
                           "Berryman1979_forced" = {
                             ## same as Berryman1979_fit
                             data.frame(
                               year = c(1:15),
                               log10Xtm1 = c(-3.1, -2.75, -2.7, -2.4, -2.3, -1.2, -1, 0.2, 0.9, 0.65,
                                             1.05, 0.95, 1.1, 1.5, 1.85),
                               log10Rt = c(0.35, 0.4, 0.1, -0.4, -0.65, 0.3, 1, 0.75, 1.2, -0.7,
                                           -0.4, 0.2, 0.45, 0.3, -0.78),
                               study = c(rep("Tunnock 1970", 9), rep("Parker 1973", 6)),
                               stringsAsFactors = TRUE
                             )
                           },
                           "Boone2011" = {
                             data <- read.csv(file.path(dataPath(sim), "BooneCurveData2.csv"))
                             data$Site <- c(rep("A", 6), rep("B", 6), rep("D", 5), rep("E", 4), rep("F", 4), rep("G", 3))
                             data$Year <- c(2000:2005, 2000:2005, 2001:2005, 2002:2005, 2002:2005, 2003:2005)
                             data
                           }
  )

  ## define growth function (from regression) for each dataset
  # mod$growthFunction <- switch(P(sim)$dataset,
  #    "Berryman1979_fit" = {
  #      function(x, s) {
  #        # TODO: check this works
  #        m <- lm(log10Rt ~ poly(log10Xtm1, 3, raw = TRUE), data = mod$growthData)
  #        s * unname(predict(m, newdata = data.frame(log10Xtm1 = x)))
  #      }
  #    },
  #    "Berryman1979_forced" = {
  #      function(x, s) {
  #        # TODO: check this works
  #        poly3.params <- c(1.1, -0.2, -0.9, -0.24)
  #        s * (poly3.params[4] * x^3 + poly3.params[3] * x^2 + poly3.params[2] * x + poly3.params[1])
  #      }
  #    },
  #    "Boone2011" = {
  #      function(x, s) {
  #        ## x is number of attacked trees (ATKTREES)
  #        ## s is scaling parameter (0,1), based on climate (CLIMATE)
  #
  #        ## mortality from emigration/dispersal
  #        # r: relative stocking value (0,1)
  #        # d: slope parameter [1,Inf)
  #        # s: scaling parameter (0,1)
  #        m_e <- function(r, d, s) {
  #          s * exp(1 - d * r)
  #        }
  #
  #        # use 2004 data as baseline for unweakened hosts (i.e., a good year for trees)
  #        m <- lm(amc::logit(PropKilled) ~ log(Attacked), data = subset(mod$growthData, Year == "2004"))
  #        a <- 0.85            ## a: slope parameter, how quickly the curve drops off
  #        d <- 3               ## d: slope parameter [1,Inf)
  #        r <- 0.2             ## r: relative stocking value (0,1) ## TODO: link this to stand depletion
  #                             ## s: scaling parameter (0,1) -- provided by climate suitability map
  #        yint2 <- 0.9         ## from MacQuarrie 2011 (Fig 3d); TODO: extract from raw data
  #        yint <- yint2 + 0.3  ## somewhat arbitrary; chosen so that the resulting curve passes 1 when flexed
  #
  #        # resulting function
  #        log(amc::hill(m$coefficients[[1]], m$coefficients[[2]], exp(a * x))) +
  #          (yint - m_e(r, d, s) - 0.03 * exp(a * x))
  #        }
  #    }
  # )
  # ! ----- STOP EDITING ----- ! #


  return(invisible(sim))
}

### plotting
plotFn <- function(sim) {
  r <- raster(sim$massAttacksStack[[1]]) ## use a template
  atkMap <- amc::dt2raster(sim$massAttacksDT, r, "ATKTREES")
  setColors(atkMap) <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
  Plot(atkMap, title = "Simulated Attacks")

  return(invisible(sim))
}

### spread
dispersal2 <- function(quotedSpread, pineMap, studyArea, massAttacksDT, massAttacksStack,
                       rasterToMatch, # massAttacksRas,
                       params, currentTime, bgSettlingProp, type, reqdPkgs,
                       windDirStack, windSpeedStack) {

  ## check that MPB and pine rasters are the same resolution and ncells
  # if (fromDisk(pineMap))
  #   pineMap[] <- pineMap[]
  #
  # ## use 1125 trees/ha, per Whitehead & Russo (2005), Cooke & Carroll (2017)
  # MAXTREES <- 1125 * prod(res(pineMap)) / 100^2 ## TODO: round this?
  #
  # ## asymmetric spread (biased eastward)
  # # lodgepole pine and jack pine together
  propPineRas <- if (nlayers(pineMap) > 1) sum(pineMap) else pineMap
  mv <- maxValue(propPineRas)
  #
  if (mv > 1)
    propPineRas[] <- propPineRas[] / 100 ## much faster than calc; ## TODO: allow different params per species

  nas <- is.na(propPineRas[]) & !is.na(rasterToMatch[])
  propPineRas[nas] <- 0

  # Put objects in Global for objFun -- this is only when not using multi-machine cluster
  advectionDir <- params$advectionDir
  advectionMag <- params$advectionMag
  omitPastPines <- FALSE # was TRUE ... but maybe bad because there are many that have multiple years
  sdDist <- 1.2
  dispersalKernel <- "Weibull"
  # dispersalKernel <- "Exponential"
  p <- do.call(c, params[c("meanDist", "advectionMag")])#, "advectionDir")])
  p["meanDist"] <- 1e4 # Average distance per year -- in m
  p["advectionMag"] <- 5e3 # Average distance of the asymmetry -- in m
  #p["advectionDir"] <- 90
  p["sdDist"] <- 1.2 # sqrt(variance(of the Weibull distribution)) --> contributes to shape and scale parameters
  #p["meanDistSD"] <- 20 # Interannual variation # not needed once we have wind and climate suitability
  #p["advectionDirSD"] <- 20
  # lower <- c(25000,  50,  0.9, 1.001)
  # upper <- c(65000, 2000, 2.5, 1.7)
  lower <- c(25000,  50,  0.9)
  upper <- c(65000, 2000, 2.5)
  p[] <- sapply(seq_along(p), function(x) runif(1, lower[x], upper[x]))

  maxDistance <- 1e5
  if (is.null(quotedSpread))
    quotedSpread <- quotedSpreadDefault
  fitType <- "distanceFromEachPoint"
  objsToExport <- setdiff(formalArgs("objFun"),
                          c("p", "reps", "quotedSpread", "fitType", "distanceFunction"))
  libPaths <- .libPaths()
  objsToExport <- c("reqdPkgs", objsToExport, "objsToExport", "libPaths")

  # RESAMPLE MAPS TO COARSER
  # nams <- names(massAttacksStack)
  # rasTmplate <- raster(massAttacksStack)
  # rasCoarse <- raster(raster::aggregate(rasTmplate, fact = 4))
  # massAttacksList <- raster::unstack(massAttacksStack)
  # opts <- options(reproducible.cacheSpeed = "fast")
  # on.exit(options(opts))
  # mams <- raster::stack(
  #   Cache(lapply, massAttacksList, function(mam)
  #     aggregateRasByDT(mam, rasCoarse, fn = sum))
  # )
  # names(mams) <- names(massAttacksStack)
  # massAttacksStack <- mams
  # propPineRas <- Cache(aggregateRasByDT, propPineRas, rasCoarse, fn = mean)
  # windDirStack <- Cache(aggregate, windDirStack, res(rasCoarse)[1]/res(windDirStack)[1])
  # windSpeedStack <- Cache(aggregate, windSpeedStack, res(rasCoarse)[1]/res(windSpeedStack)[1])
  # options(opts)

  list2env(mget(objsToExport), envir = .GlobalEnv)

  ips <- c("localhost", "10.20.0.184", #"10.20.0.97",
           "10.20.0.220", "10.20.0.217") ## TODO: don't hardcode this. pass as param
  #adjustments <- c(1.25, 1.1, #0.8,
  #                 0.8, 1.4) # relative, manual, highly idiosyncratic, depends on current use of machines
  #ips <- c("10.20.0.106", "10.20.0.68", "10.20.0.213")
  adjustments = c(0.9, 1.2, 0.8, 1.2)
  if (isTRUE(type == "fit")) {
    fn <- ".allObjs.rda"
    numCoresNeeded <- 84 # This is one per true core on 4 machines, 56, 56, 28, 28 cores each
    reqdPkgs <- grep(paste(collapse = "|", c("SpaDES.tools", "raster", "CircStats", "data.table", "purrr")), reqdPkgs, value = TRUE)
    clusterIPs <- clusterSetup(workers = ips, objsToExport = objsToExport,
                               packages = reqdPkgs, libPaths = libPaths, doSpeedTest = TRUE, fn = fn,
                               numCoresNeeded = numCoresNeeded, adjustments = adjustments)
    cl <- clusterSetup2(clusterIPs, filenameToLoad = fn, packages = reqdPkgs)
    on.exit(try(stopCluster(cl)), add = TRUE)
    message("Starting DEoptim")
    stPre <- Sys.time()
    DEout <- DEoptim(fn = objFun,
                     lower = lower,#c(500, 1, -90, 0.9, 1.1, 5),
                     upper = upper,#c(30000, 10, 240, 1.8, 1.6, 30),
                     reps = 1,
                     quotedSpread = quotedSpread,
                     control = DEoptim.control(cluster = cl,
                                               strategy = 6,
                                               itermax = 200,
                                               NP = numCoresNeeded),
                     fitType = fitType)

    saveRDS(DEout, file = file.path("outputs", paste0("DEout_", format(stPre), ".rds")))
    sim2 <- get("sim", whereInStack("sim"))
    parms <- params(sim2)
    saveRDS(parms, file = file.path("outputs", paste0("parms_", format(stPre), ".rds")))
    stPost <- Sys.time()
    print(difftime(stPost, stPre))
    try(parallel::stopCluster(cl))
    browser()

  } else if (isTRUE(type == "runOnce")) {
    #profvis::profvis(
      st1 <- system.time({
        out <- objFun(quotedSpread = quotedSpread, reps = 1, p = p, fitType = fitType,
                      massAttacksStack = massAttacksStack, massAttacksDT = massAttacksDT)
      })
    #)
    print(st1)
  } else if (isTRUE(type == "optim")) {
    browser()
    cl <- parallel::makeForkCluster(5)
    stPre <- Sys.time()
    optimOut <- optimParallel::optimParallel(fn = objFun, par = p,
                     lower = lower,#c(500, 1, -90, 0.9, 1.1, 5),
                     upper = upper,#c(30000, 10, 240, 1.8, 1.6, 30),
                     reps = 1,
                     method = "L-BFGS-B",
                     quotedSpread = quotedSpread,
                     control = list(trace = 3, factr = 1e6),
                     parallel = list(cl = cl, forward = TRUE, loginfo  = TRUE),
                     #control = DEoptim.control(cluster = cl,
                    #                           strategy = 6,
                    #                           itermax = 120,
                    #                           NP = numCoresNeeded),
                     fitType = fitType)
    saveRDS(optimOut, file = file.path("outputs", paste0("optimOut_", format(stPre), ".rds")))
    stPost <- Sys.time()

  }
  browser()
  objs <- dir("outputs", pattern = "DEout", full.names = TRUE)
  DEout <- readRDS(tail(objs, 1))
  p[] <- DEout$optim$bestmem[] # use [] to keep names


  N1 <- 1e4

  plotKernels(p)

  #########################################################
  #################### PREDICT MAPS
  #########################################################

  opts2 <- options(reproducible.cacheSpeed = "fast")
  on.exit(options(opts2))
  # Browse[1]> fastdigest(quotedSpread)
  # [1] "6ee73bef9a187c2404d856d22b7a426a"
  predictedDT <- Cache(predictQuotedSpread,
                       massAttacksStack = massAttacksStack,
                       windDirStack = windDirStack,
                       windSpeedStack = windSpeedStack,
                       propPineRas = propPineRas,
                       dispersalKernel = "Weibull",
                       maxDistance = maxDistance,
                       quotedSpread = quotedSpread, # doesn't cache correctly
                       quotedSpreadChar = format(quotedSpread), # cache this instead
                       p = p,
                       omitArgs = "quotedSpread"
  )

  predictedStack <- Cache(stackFromDT, massAttacksStack, predictedDT)

  plotStack <- Cache(stacksForPlot, massAttacksStack, predictedStack, threshold = 0.5)

  browser()
  return(invisible(plotStack))
  browser()
  #########################################################
  #########################################################

  visualizeKernel(p, massAttacksStack, maxDistance = maxDistance)
  ###############################################

  migrantsDT <- out[, list(NEWATKs = sum(abundSettled)), by = "pixels"]
  massAttacksDTNew <- massAttacksDT[migrantsDT, on = c(ID = "pixels")]

  # ! ----- STOP EDITING ----- ! #
  return(invisible(massAttacksDTNew))
}



spiralDistances <- function(pixelGroupMap, maxDis, cellSize) {
  spiral <- which(focalWeight(pixelGroupMap, maxDis, type = "circle")>0, arr.ind = TRUE) -
    ceiling(maxDis/cellSize) - 1
  spiral <- cbind(spiral, dists = sqrt( (0 - spiral[,1]) ^ 2 + (0 - spiral[, 2]) ^ 2))
  spiral <- spiral[order(spiral[, "dists"], apply(abs(spiral), 1, sum),
                         abs(spiral[, 1]), abs(spiral[, 2])),, drop = FALSE]
}

distanceFunction <- function(dist, angle, landscape, fromCell, toCells, # nextYrVec,
                             propPineRas, asymParam, windDir, meanDist, sdDist,
                             dispersalKernel, maxDistance, windSpeed) {

  if (asymParam > 0) {
    fromXY <- xyFromCell(propPineRas, fromCell)
    toXYs <-  xyFromCell(propPineRas, toCells)
    if (length(windDir) > 1) {
      windDir <- windDir[fromCell]
      windSpeed <- windSpeed[fromCell]
    }
    advectionXY <- c(x = sin(rad(windDir))*asymParam*windSpeed, y = cos(rad(windDir))*asymParam*windSpeed)
    x <- toXYs[, 'x'] - advectionXY[1]/meanDist*dist
    y <- toXYs[, 'y'] - advectionXY[2]/meanDist*dist
    newTo <- cbind(x, y)
    colnames(newTo) <- c("x", "y")
    if (!is.null(dim(dist))) dist <- dist[,1, drop = TRUE]
    dists <- cbind(origX = toXYs[, "x"], origY = toXYs[, "y"],
                   origDist = dist,
                   .pointDistance(from = fromXY, newTo, angles = FALSE))
    dist <- dists[, "dists"]
  }

  prob <- kernelFn(meanDist, sdDist, dispersalKernel, dist)

  ll <- landscape[fromCell]
  prob1 <- prob * ll
  return(prob1)
}

aggregateRasByDT <- function(ras, newRas, fn = sum) {
  whNonNA <- which(!is.na(ras[]))
  rc2 <- rowColFromCell(ras, whNonNA)
  if (!all(((res(newRas)/res(ras)) %% 1) == 0))
    stop("The resolutions of the original raster and new raster are not integer multiples")
  disaggregateFactor <- unique(res(newRas)/res(ras))
  dt <- data.table(vals = ras[][whNonNA], ceiling(rc2 / disaggregateFactor))
  dt2 <- dt[, list(vals = fn(vals)), by = c("row", "col")]
  pixes <- cellFromRowCol(newRas, row = dt2$row, col = dt2$col)
  newRasOut <- raster(newRas)
  newRasOut[pixes] <- dt2$vals
  newRasOut

}

clusterSetup <- function(workers, objsToExport, packages,
                         libPaths, doSpeedTest = FALSE, envir = parent.frame(),
                         fn = ".allObjs.rda", numCoresNeeded, adjustments = rep(1, length(workers))) {
  message("Starting cluster with 1 core per machine -- install packages; copy objects; write to disk")
  clSingle <- future::makeClusterPSOCK(workers = workers, revtunnel = TRUE)
  on.exit(try(parallel::stopCluster(clSingle), silent = TRUE), add = TRUE)
  # NumPopulations <- 118

  if (isTRUE(doSpeedTest)) {
    message("testing speed on each to estimate number cores to use")
    out <- clusterEvalQ(clSingle, {
      ss <- system.time({
        for (i in 1:10000) rnorm(1e4)
      })
      data.table::data.table(elapsed = ss[3], trueCores = parallel::detectCores())
    })
  } else {
    out <- clusterEvalQ(clSingle, {
      data.table::data.table(elapsed = 1, trueCores = parallel::detectCores())
    })
  }

  # parallel::clusterEvalQ(clSingle, system("pkill -f workRSOCK"))
  names(out) <- workers
  out2 <- rbindlist(out, idcol = TRUE)
  relSpeed <- out2$trueCores/out2$elapsed*numCoresNeeded
  relSpeed <- relSpeed/numCoresNeeded
  relSpeed <- relSpeed/max(relSpeed)
  out2[, relSpeed := relSpeed]
  nonHTcores <- out2$trueCores/2
  out2[, nonHTcores := nonHTcores]
  sumNonHTcores <- sum(nonHTcores)
  # needHTcores <- max(numCoresNeeded, numCoresNeeded - sumNonHTcores)

  out2[, cores := round(nonHTcores * adjustments / max(adjustments))]

  # Increase ncores upwards
  m <- 0
  while (sum(out2$cores) < numCoresNeeded) {
    m <- ((m + 1) - 1) %% NROW(out2) + 1
    out2[m, cores := cores + relSpeed]
  }
  # Decrease down to 1 each
  out2[, cores := ceiling(cores)]
  m <- 0
  set(out2, NULL, "cores", as.numeric(out2$cores))
  while ((sum(floor(out2$cores)) > numCoresNeeded) && any(out2$cores > 1)) {
    m <- ((m + 1) - 1) %% NROW(out2) + 1
    if (out2[m,]$cores > 1) {
      maxC <- max(out2$cores)
      out2[m, cores := cores - 1]
    }
  }
  # set(out2, NULL, "cores", floor(out2$cores))
  # Decrease down to 0 or 1 each
  m <- 0
  while (sum(out2$cores) > numCoresNeeded && any(out2$cores > 0)) {
    m <- ((m + 1) - 1) %% NROW(out2) + 1
    if (out2[m,]$cores > 0) {
      out2[m, cores := cores - 1]
    }
  }

  reproducible::messageDF(out2)

  clusterIPs <- rep(out2$.id, out2$cores)
  # parallel::stopCluster(clSingle)

  # clSingle <- future::makeClusterPSOCK(workers = unique(clusterIPs), revtunnel = TRUE)
  clusterExport(clSingle, varlist = objsToExport, envir = envir)
  clusterExport(clSingle, varlist = c("fn", "packages", "libPaths", "objsToExport"), envir = environment())
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

    suppressMessages(Require::Require(packages, install = TRUE, require = FALSE))
    save(list = objsToExport, file = fn)
  })
  parallel::stopCluster(clSingle)
  return(clusterIPs)
  #clOut <- future::makeClusterPSOCK(workers = clusterIPs, revtunnel = TRUE)
  #return(clOut)
}


grow <- function(sim) {
  ## determine the actual growth based on the actual number of attacked trees/ha
  sim$massAttacksDT[, greenTreesYr_t := xt(ATKTREES, CLIMATE, P(sim)$dataset,
                                           sim$massAttacksStack, mod$growthData)]

  return(invisible(sim))
}

# clOut= clusterSetup(workers = ips[1], numCoresNeeded = 1,
#                     objsToExport = c("ips"), packages = "reproducible", libPaths= .libPaths()[1])

hist.DEoptim <- function(DEobj, paramNames) {
  dt <- as.data.table(DEobj$member$pop)
  colnames(dt) <- paramNames
  dt <- melt(dt, measure.vars = seq(NCOL(dt)))
  gg <- ggplot(dt, aes(x = value)) +
    geom_histogram(bins = 15) +
    facet_grid(. ~ variable, scales = "free") +
    theme_bw()

  gg
}

summary.DEoptim <- function(DEobj, title) {
  if (!missing(title))
    cat(paste(title, "\n"))
  cat(paste("Likelihood: ", round(DEobj$optim$bestval, 0), "\n"))
  cat(paste("Mean params: ", paste(round(apply(DEobj$member$pop, 2, mean), 2), collapse = ", "), "\n"))
  cat(paste("Best params: ", paste(round(DEobj$optim$bestmem, 2), collapse = ", "), "\n"))
}

clusterSetup2 <- function(ips, packages, filenameToLoad) {
  message("Starting cluster with all cores per machine")
  cl <- future::makeClusterPSOCK(ips, revtunnel = TRUE)
  # cl <- future::makeClusterPSOCK(workers = clusterIPs, revtunnel = TRUE)
  # on.exit(try(parallel::stopCluster(cl)), add = TRUE)
  clusterExport(cl, varlist = c("filenameToLoad", "packages"), envir = environment())
  st <- system.time(clusterEvalQ(cl, {
    load(file = filenameToLoad, envir = .GlobalEnv)
    .libPaths(libPaths)
    suppressMessages(
      Require::Require(packages, install = FALSE)
    )
  }))
  st <- system.time(clusterEvalQ(cl, {
    try(unlink(filenameToLoad), silent = TRUE)
  }))

  return(cl)
}

plotKernels <- function(p, reps = 10) {
  nrows <- max(1, floor(sqrt(reps + 1)))
  ncols <- max(1, ceiling((reps + 1)/nrows))
  par(mfrow = c(nrows, ncols))
  if (is.null(names(p))) stop("p must be named and have these 3: meanDist, sdDist, meanDistSD")
  for (i in 1:reps) {
    meanDist <- rlnorm(1, log(p["meanDist"]), log(p["meanDistSD"]))
    sdDist = p["sdDist"]

    mn <- (meanDist)
    sd <- mn/sdDist # 0.8 to 2.0 range
    shape <- (sd/mn)^(-1.086)
    scale <- mn/exp(lgamma(1 + 1/shape))

    ww <- rweibull(1e4, shape = shape, scale = scale)
    hist(ww, xlim = c(0, 14e4), main = "", xlab = "Distance (m) from source")
    abline(v = meanDist, col = "red")
    abline(v = p[[1]], col = "green")
  }
  plot(type = "n", x = 1:10, axes = F, xlab = "", ylab = "")
  legend(x = 2, y = 8, col = c("red", "green"), lty= 1,
         legend = c("annual average", "global average"),
         cex = 1.5)
  mtext(text = "Dispersal kernels", outer = TRUE, line = -3, xpd = TRUE)
}

kernelFn <- function(meanDist, sdDist, dispersalKernel, dist) {
  if (grepl("Weibull", dispersalKernel)) {
    mn <- (meanDist)
    sd <- mn/sdDist # 0.8 to 2.0 range
    shape <- (sd/mn)^(-1.086)
    scale <- mn/exp(lgamma(1+1/shape))

    if (grepl("Weibull", dispersalKernel)) {
      prob <- dweibull(dist, scale = scale, shape = shape)
    }

    infs <- is.infinite(prob)
    if (any(infs))
      prob[infs] <- 0
  } else {
    prob <- dexp(dist, rate = 1/meanDist)

  }
  prob
}


stackFromDT <- function(massAttacksStack, predictedDT) {
  indices <- matrix(rep(NA, nlayers(massAttacksStack) * ncell(massAttacksStack)),
                    ncol = nlayers(massAttacksStack))
  yrNames <- names(massAttacksStack)
  yrNamesPlus1 <- yrNamesPlus1(yrNames)
  colnames(indices) <- yrNamesPlus1


  #saveFilename <- file.path("outputs", paste0("distanceSurface_", nams[i], "_", Sys.time(), ".qs"))
  #qs::qsave(predictedDT, file = saveFilename)
  pixels <- cellFromXY(massAttacksStack, predictedDT[, c("x", "y")])
  predictedStack <- massAttacksStack
  set(predictedDT, NULL, "cell", cellFromXY(massAttacksStack, predictedDT[, c("x", "y")]))
  outW <- dcast(predictedDT[, c("Year", "cell", "val")], cell  ~ Year)
  cols <- grep("^X", colnames(outW), value = TRUE)
  ff <- as.integer(factor(predictedDT$Year))
  indices[cbind(predictedDT$cell, ff)] <- as.matrix(outW[, ..cols])

  predictedStack[] <- indices
  predictedStack <- raster::stack(predictedStack)
  predictedStack
}

yrNamesPlus1 <- function(yrNames) {
  yrNames <- gsub("X([[:digit:]]{4,4})", "\\1", yrNames)
  paste0("X", as.integer(yrNames) + 1)
}

stacksForPlot <- function(massAttacksStack, predictedStack, threshold = 1) {
  yrNames <- names(massAttacksStack)
  yrNamesPlus1 <- yrNamesPlus1(yrNames)
  stk <- raster::stack()
  for (lay in seq(yrNames)) {
    setColors(predictedStack[[yrNamesPlus1[lay]]], 8) <- "Reds"
    AllThree <- raster(massAttacksStack[[lay]])
    whCA <- which(!is.na(massAttacksStack[[lay]][]))
    AllThree[] <- 0
    AllThree[whCA] <- 1
    if (yrNamesPlus1[lay] %in% names(massAttacksStack)) {
      whNY <- which(!is.na(massAttacksStack[[yrNamesPlus1[lay]]][]))
      AllThree[whNY] <- AllThree[whNY] + 2
    }
    whPM <- which(!is.na(predictedStack[[yrNamesPlus1[lay]]][]) &
                    predictedStack[[yrNamesPlus1[lay]]][] > threshold)
    AllThree[whPM] <- AllThree[whPM] + 4
    levels(AllThree) <- data.frame(ID = 1:7,
                                   Label = c("LastYr", "NextYr", "LastAndNext", "PredYr", "LastYrAndPredYr", "NextYrAndPredYr", "LYRandNYandPY"))
    setColors(AllThree, n = 7) <- "Set2"
    AllThree@legend@colortable[6:7] <- "red"
    stk <- addLayer(stk, AllThree)
  }
  names(stk) <- paste0("pred", yrNamesPlus1)
  stk[which(stk[] == 0)] <- NA
  stk
}

