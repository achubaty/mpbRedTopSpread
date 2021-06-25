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
                  "gamlss", "ggplot2", "gt",
                  "PredictiveEcology/SpaDES.core@modifyList2 (>= 1.0.8.9011)",
                  "PredictiveEcology/LandR@cluster (>= 1.0.4.9022)",
                  "parallelly",
                  "PredictiveEcology/pemisc@development (>= 0.0.3.9001)",
                  "PredictiveEcology/mpbutils (>= 0.1.3)", "purrr",
                  "quickPlot", "raster", "RColorBrewer",
                  "PredictiveEcology/reproducible@DotsBugFix (>= 1.2.7.9009)",
                  "PredictiveEcology/SpaDES.tools@development (>= 0.3.7.9021)",
                  "tmap"),
  parameters = rbind(
    defineParameter("cachePredict", "logical", TRUE, NA, NA,
                    "The function predictQuotedSpread can be Cached or not; default is TRUE"),
    defineParameter("dataset", "character", "Boone2011", NA, NA, "Which dataset to use for stand dynamic model fitting. One of 'Boone2011' (default), 'Berryman1979_fit', or 'Berryman1979_forced'. Others to be implemented later."),
    defineParameter("growthInterval", "numeric", 1, NA, NA, "This describes the interval time between growth events"),
    # defineParameter("p_advectionDir", "numeric", 90, 0, 359.9999,
    #                 "The direction of the spread bias, in degrees from north"),
    # defineParameter("bgSettlingProp", "numeric", 0.1, 0, 1,
    #                 "The proportion of beetles that settle from those that could potentially settle, even if no pine"),
    defineParameter("dispersalKernel", "character", "Weibull3", "GeneralizedGamma", "Exponential",
                    "This can also be a character string: one of Weibull3, GeneralizedGamma, twodt, or Exponential to use one of these
                    different dispersal kernels."),
    defineParameter("maxDistance", "numeric", 1.4e5, NA, NA,
                    "The maximum distance to allow for pair-wise from-to dispersal pairs"),
    defineParameter("quotedSpread", c("language"), NULL, NA, NA,
                    "A spread function that is contained within a quote(). If left at NULL or NA, it will use the module version."),
    defineParameter("p_advectionMag", "numeric", 5e3, NA, NA,
                    "The magnitude of the wind effect on spread. This number is multiplied by the wind speed (which averages 9.6 km/h in the historical dataset)"),
    defineParameter("p_meanDist", "numeric", 1e4, NA, NA,
                    "Expected dispersal distance (m); ~63% go less than this distance"),
    defineParameter("p_nu", "numeric", 1, NA, NA,
                    "Third parameter in the generalized gamma"),
    # defineParameter("p_rescaler", "numeric", 1, NA, NA,
    #                 "Single scalar that rescales the predicted abundances that result from growFit ",
    #                 "& dispersalFit, so that they match the observed data"),
    defineParameter("p_sdDist", "numeric", 1.2, NA, NA,
                    paste0("The dispersion term for the Weibull dispersal kernel: contributes to shape and scale parameters;",
                           "sqrt(variance(of the Weibull distribution)) ")),
    defineParameter("type", "character", "DEoptim", NA, NA,
                    "One of several modes of running this module: DEoptim, optim, runOnce, validate, predict or nofit"),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Zero or more of c('screen', 'png', 'pdf', 'raw'. See ?Plots"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
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
    expectsInput(objectName = "absk", objectClass = "SpatialPolygonsDataFrame",
                 desc = "Alberta and Saskatchewan political outlines"),
    expectsInput("climateSuitabilityMaps", "RasterStack",
                 "A time series of climatic suitablity RasterLayers, each with previx 'X' and the year, e.g., 'X2010'"),
    expectsInput("massAttacksStack", "RasterStack",
                 desc = "Historical MPB attack maps (number of red attacked trees).",
                 sourceURL = NA),
    # expectsInput("massAttacksDT", "data.table",
    #              desc = "Current MPB attack map (number of red attacked trees).",
    #              sourceURL = NA),
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
    createsOutput("massAttacksStack", "RasterStack",
                 "This will be the same data as the inputted object, but will have a different resolution."),
    createsOutput("fit_mpbSpreadOptimizer", "list", "The output object from DEoptim"),
    createsOutput("propPineRas", "RasterLayer",
                  paste("Proportion (not percent) cover map of *all* pine. This will",
                        " be the pixel-based sum if there are more than one layer")),
    createsOutput("predictedDT", "data.table",
                  "Similar to massAttacksDT, but with predicted (from the models here) values of ATKTREES"),
    createsOutput("thresholdAttackTreesMinDetectable", "numeric",
                  paste("The number of predicted attacks (post dispersal) per pixel below which ",
                        "it can be considered 'no attack'")),
    createsOutput("ROCList", "list",
                  "")


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
           if (isTRUE(any(c("fit", "DEoptim", "runOnce", "optim", "all") %in% P(sim)$type))) {
             sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "growFit", eventPriority = 2.5)
             sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "dispersalFit", eventPriority = 3.5)
           }
           if (isTRUE(any(c("validate", "all") %in% P(sim)$type))) {
             sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "growFit", eventPriority = 2.5)
             sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "validate", 4.5)
           }
           if (isTRUE(any(c("predict", "all") %in% P(sim)$type))) { # will automatically validate
             sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "growPredict", eventPriority = 5.5)
             sim <- scheduleEvent(sim, end(sim), "mpbRedTopSpread", "validate", 7)
           } else if (isTRUE(any(c("validateAll") %in% P(sim)$type))) {
             sim <- scheduleEvent(sim, end(sim), "mpbRedTopSpread", "validate", 7)
           }

           # sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "plot", 7.5)
           # sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "mpbRedTopSpread", "save", .last())
         },
         "growFit" = {
           sim$massAttacksDT <- grow(sim$massAttacksDT, P(sim)$dataset,
                                     sim$massAttacksStack, mod$growthData)
         },
         "growPredict" = {
           if (is.null(sim$fit_mpbSpreadOptimizer)) {
             message("A fit_mpbSpreadOptimizer object should probably be supplied")
             optimOutFileList <- dir("outputs", pattern = "optimOut", full.names = TRUE)
             sim$fit_mpbSpreadOptimizer <- readRDS(optimOutFileList[[3]])
             DEoutFileList <- dir("outputs", pattern = "fit_mpbSpreadOptimizer", full.names = TRUE)
             sim$fit_mpbSpreadOptimizer <- readRDS(DEoutFileList[[23]])
             sim$fit_mpbSpreadOptimizer <- colnamesToDEout(sim$fit_mpbSpreadOptimizer, c("p_meanDist", "p_advectionMag", "p_sdDist"))
           }
           sim$massAttacksDT_1Yr <- growPredict(sim)

           ## growth needs to happen after spread:
           sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "dispersalPredict", eventPriority = 5.5)
         },
         "dispersalPredict" = {
           sim$predictedDT <- dispersalPredict(sim) # this object accumulates years over time
           sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "growPredict", eventPriority = 5.5)
         },
         "dispersalFit" = {
           sim$fit_mpbSpreadOptimizer <- dispersalFit(quotedSpread = P(sim)$quotedSpread,
                             propPineRas = sim$propPineRas, studyArea = sim$studyArea,
                             massAttacksDT = sim$massAttacksDT,
                             massAttacksStack = sim$massAttacksStack,
                             maxDistance = P(sim)$maxDistance,
                             params = P(sim),
                             windDirStack = sim$windDirStack,
                             windSpeedStack = sim$windSpeedStack,
                             rasterToMatch = sim$rasterToMatch,
                             # bgSettlingProp = P(sim)$bgSettlingProp,
                             dispersalKernel = P(sim)$dispersalKernel,
                             type = P(sim)$type,
                             # currentTime = time(sim),
                             reqdPkgs = reqdPkgs(module = currentModule(sim), modulePath = modulePath(sim))[[currentModule(sim)]]
           )

           # sim <- dispersal(sim)
           # sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "dispersal", eventPriority = 5.5)
         },
         "validate" = {
           sim <- Validate(sim)
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

  kern <- tolower(P(sim)$dispersalKernel)
  if (grepl("weibul", kern))
    P(sim, "quotedSpread") <- quotedSpreadWeibull3
  if (grepl("exponen", kern))
    P(sim, "quotedSpread") <- quotedSpreadExponential
  if (grepl("gamma", kern))
    P(sim, "quotedSpread") <- quotedSpreadGeneralGamma

  if (is.null(P(sim)$quotedSpread))
    P(sim, "quotedSpread") <- quotedSpreadWeibull3


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
  mv <- maxValue(sim$pineMap)
  propPineRas <- raster(sim$pineMap)
  propPineRas[] <- if (mv > 1)
    sim$pineMap[] <- sim$pineMap[] / 100 ## much faster than calc; ## TODO: allow different params per species
  else
    sim$pineMap[]

  nas <- is.na(propPineRas[])# & !is.na(rasterToMatch[])
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
  ),
  by = "layerName"]

  sim$predictedDT <- sim$massAttacksDT[0]

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

  return(invisible(sim))
}

dispersalPredict <- function(sim) {

  pBest <- sim$fit_mpbSpreadOptimizer$optim$bestmem

  colNameForPrediction <- "ATKTREES"
  predictedDT <- Cache(predictQuotedSpread,
                       massAttacksDT = sim$massAttacksDT_1Yr,
                       windDirStack = sim$windDirStack,
                       windSpeedStack = sim$windSpeedStack,
                       propPineRas = sim$propPineRas,
                       pineThreshold = 0.3,
                       dispersalKernel = P(sim)$dispersalKernel,
                       clNumber = 12,
                       useCache = P(sim)$cachePredict,
                       maxDistance = P(sim)$maxDistance,
                       quotedSpread = P(sim)$quotedSpread, # doesn't cache correctly
                       .cacheExtra = format(P(sim)$quotedSpread), # cache this instead
                       p = pBest,
                       # colNameForPrediction = colNameForPrediction,
                       omitArgs = "quotedSpread"
  )

  # predictedDT <- predictedDT[get(colNameForPrediction) > sim$thresholdAttackTreesMinDetectable]
  predictedDT[, CLIMATE := sim$climateSuitabilityMaps[[unique(layerName)]][pixel]]
  return(rbindlist(list(sim$predictedDT, predictedDT), use.names = TRUE, fill = TRUE))
}

Validate <- function(sim) {
  startToEnd <- paste0(start(sim), " to ", end(sim), "_", Sys.time())

  if (is.null(sim$DEout)) {
    DEoutFileList <- dir("outputs", pattern = "DEout", full.names = TRUE)
    fit_mpbSpreadOptimizer <- readRDS(tail(DEoutFileList, 1)) # 24 was best so far
  }

  DEouts <- list(fit_mpbSpreadOptimizer)

  sim$ROCList <- list()

  for (iii in seq_along(DEouts)) {#seq(DEoutFileList)) { # 22 is good (the 6 parameter generalgamma)
    if (FALSE) {
      optimOut <- readRDS(optimOutFileList[3])
      p <- optimOut$par

      objsToExport <- setdiff(formalArgs("objFun"),
                              # The following are small and passed at the time of the DEoptim call
                              c("p", "quotedSpread", "distanceFunction"))
      sim$omitPastPines <- FALSE
      sim$p_sdDist <- p["p_sdDist"]
      sim$dispersalKernel <- P(sim)$dispersalKernel
      sim$maxDistance <- Par$maxDistance
      ss <- do.call(numDeriv::hessian, append(list(func = objFun),
                                              mget(objsToExport, envir = envir(sim))))

    }

    if (length(fit_mpbSpreadOptimizer$optim$bestmem) %in% 3:4 || iii > 20) {
      if (all(grepl("par[[:digit:]]", names(fit_mpbSpreadOptimizer$optim$bestmem))))
        if (length(fit_mpbSpreadOptimizer$optim$bestmem) == 3) {
          fit_mpbSpreadOptimizer <- colnamesToDEout(fit_mpbSpreadOptimizer, c("p_meanDist", "p_advectionMag", "p_sdDist"))
        } else if  (length(fit_mpbSpreadOptimizer$optim$bestmem) == 4) {
          fit_mpbSpreadOptimizer <- colnamesToDEout(fit_mpbSpreadOptimizer, c("p_meanDist", "p_advectionMag", "p_sdDist", "p_rescaler"))
        }

      p <- fit_mpbSpreadOptimizer$optim$bestmem # replace with values from disk object


      DEoutPop <- as.data.table(fit_mpbSpreadOptimizer$member$pop)
      setnames(DEoutPop, new = names(p))
      DEoutPop <- melt(DEoutPop, measure = colnames(DEoutPop), variable.name = "Parameter")

      Plots(DEoutPop, plot_HistsOfDEoptimPars, title = "Parameter histograms from DEoptim",
            filename = paste0("Histograms of parameters from DEoptim ", Sys.time()))


      # Way to identify quantiles
      if (FALSE) {
        ll <- weibullShapeScale(p['p_meanDist'], p['p_sdDist']);
        qweibull(0.99, shape = ll$shape, scale = ll$scale)
      }

      #########################################################
      #################### PREDICT MAPS
      #########################################################

      opts2 <- options(reproducible.cacheSpeed = "fast")
      on.exit(options(opts2))
      if (NROW(sim$predictedDT) == 0)
        sim$predictedDT <- Cache(predictQuotedSpread,
                                 massAttacksDT = sim$massAttacksDT, # must have column named greenTreesYr_t
                                 windDirStack = sim$windDirStack,
                                 windSpeedStack = sim$windSpeedStack,
                                 propPineRas = sim$propPineRas,
                                 pineThreshold = 0.3,
                                 dispersalKernel = P(sim)$dispersalKernel,#  "generalgamma",
                                 maxDistance = P(sim)$maxDistance,
                                 quotedSpread = P(sim)$quotedSpread, # doesn't cache correctly
                                 .cacheExtra = format(P(sim)$quotedSpread), # cache this instead
                                 p = p,
                                 omitArgs = "quotedSpread"
        )
      sim$predictedStack <- Cache(stackFromDT, sim$predictedDT, raster(sim$massAttacksStack))
      #bb <- predictedStack$X2011[][propPineRas[] < 0.3]
      #cc <- massAttacksStack$X2011[][propPineRas[] < 0.1]

      # lapply(seq(nlayers(sim$massAttacksStack)), function(x)
      #   cor.test(sim$predictedStack[][, x], sim$massAttacksStack[][, x], use = "complete.obs"))

      if (length(names(sim$massAttacksStack)) > 1) {
        corByYr <- cor(sim$predictedStack[], sim$massAttacksStack[], use = "complete.obs")
        message("Avg correlation of time series: ", round(corByYrAvg <- mean(diag(corByYr)), 3))
      }

      plotStack <- Cache(createStackPredAndObs, sim$massAttacksStack, sim$predictedStack, threshold = 55)

      fnHash <- fastdigest(format(stacksPredVObs))

      availablePine <- rowSums(sim$massAttacksStack[], na.rm = TRUE)
      whAttacked <- which(availablePine > 0)
      availablePineMap <- raster(sim$massAttacksStack)
      availablePineMap[whAttacked] <- 1
      apmBuff <- Cache(raster::buffer, availablePineMap, width = 1e5)
      apmBuff[sim$propPineRas[] == 0] <- NA

      # cumulative maps
      mam <- sim$massAttacksStack
      for (lay in names(sim$massAttacksStack)) {
        mam[[lay]][is.na(mam[[lay]][])] <- 0
      }
      massAttacksStackLog <- if (nlayers(mam) > 1) sum(mam) else mam[[1]]
      massAttacksStackLog[] <- log(massAttacksStackLog[] + 1)

      # cumulative maps
      pred <- sim$predictedStack
      for (lay in names(pred)) {
        pred[[lay]][is.na(pred[[lay]][])] <- 0
      }
      predictCumulativeLog <- sum(pred)
      predictCumulativeLog[] <- log(predictCumulativeLog[] + 1)


      if (tolower(Par$type) == "validate") {
        # Find numAttackTrees minimum threshold that defines detectable
        message("Estimating (using ROC) the predicted attack density minimum threshold ",
                " to be considered detectable")
        system.time(optimAttackTreesMinDetectable <- Cache(optimize, f = estimateMinAttackThresh, interval = c(0.001, 200),
                                                           meanAttackStk = sim$massAttacksStack, .cacheExtra = fnHash,
                                                           predictedStk = sim$predictedStack, propPineRast = sim$propPineRas,
                                                           tol = 0.001
        ))
        sim$thresholdAttackTreesMinDetectable <- optimAttackTreesMinDetectable$minimum
      } else {
        sim$thresholdAttackTreesMinDetectable <- 4.0 # this is best from DEoptim
      }


      # ROC curves
      stackPredVObs <- Cache(stacksPredVObs, sim$massAttacksStack, sim$predictedStack,
                                      apmBuff, threshold = sim$thresholdAttackTreesMinDetectable,
                                      plot = TRUE)

      obsPredDT <- makeAnnualMeansDT(sim$predictedStack, sim$massAttacksStack)
      Plots(fn = plot_ObsVPredAbund1, dt = obsPredDT, #types = "screen",
            filename = paste0("Annual Means: Obs v Pred ", startToEnd))
      Plots(fn = plot_ObsVPredAbund1, dt = obsPredDT, #types = "screen",
            filename = paste0("Annual Means: Year v Obs and Pred ", startToEnd))

      # browser()
      print(paste("ROC mean value: ", round(mean(stackPredVObs$ROCs$ROC, na.rm = TRUE), 3)))

      messageDF(stackPredVObs$ROCs)
      sim$ROCList[[iii]] <- stackPredVObs$ROCs
    }
  }
  # current.mode <- tmap_mode("view")

  current.mode <- tmap_mode("plot")
  clearPlot()

  fn1 = file.path(outputPath(sim), "figures", paste0(Par$type, ": stks predicted vs. observed ", startToEnd))

  # tmap doesn't work with Plots yet... do manually
  #Plots(sim$absk, plot_StudyAreaFn,
  #filename = file.path(outputPath(sim), paste0(Par$type, ": stks predicted vs. observed ", startToEnd)),
  #ggsaveArgs = list(width = 8, height = 10, units = "in", dpi = 300),
  tm_stks_predvObs <- plot_StudyAreaFn(
        absk = sim$absk, cols = "darkgreen", sf::st_as_sf(sim$studyArea),
        massAttacksStackLog, pred = predictCumulativeLog,
        propPineMap = sim$propPineRas,
        minThreshold = sim$thresholdAttackTreesMinDetectable)
  if ("screen" %in% Par$.plots) {
    tm_stks_predvObs
  }
  if (any(Par$.plots %in% ggplotClassesCanHandle)) {
    tmap_save(tm_stks_predvObs, filename = paste0(fn1, ".", Par$.plots),
              width = 8, height = 10, units = "in", dpi = 300)
  }

  ss <- lapply(raster::unstack(stackPredVObs$stk), function(ras) {
    freqs <- table(ras[])
    sensitivity <- freqs[3]/(sum(freqs[3], freqs[1]))
    specificity <- freqs[4]/(sum(freqs[4], freqs[2]))
    return(list(sensitivity = sensitivity, specificity = specificity))
  })


  # Edges
  edgesPredicted <- sim$predictedStack
  edgesPredicted[edgesPredicted < sim$thresholdAttackTreesMinDetectable] <- NA
  edgesPredicted <- raster::stack(edgesPredicted)
  edgesList <- raster::unstack(edgesPredicted)
  names(edgesList) <- names(edgesPredicted)
  edgesCoords <- lapply(edgesList, function(edge) {
    xys <- xyFromCell(edge, which(!is.na(edge[])))
    xEdge <- quantile(xys[, "x"], 0.95)
    yEdge <- quantile(xys[, "y"], 0.95)
    list(EastEdge = xEdge, NorthEdge = yEdge)
  })
  edgesPredDT <- rbindlist(edgesCoords, idcol = "layerName")
  edgesPredDT[, type := "Prediction"]

  edgesData <- sim$massAttacksStack
  # edgesData <- raster::stack(edgesData)
  edgesDataList <- raster::unstack(edgesData)
  names(edgesDataList) <- names(edgesData)
  edgesDataCoords <- lapply(edgesDataList, function(edge) {
    xys <- xyFromCell(edge, which(!is.na(edge[])))
    xEdge <- quantile(xys[, "x"], 0.95)
    yEdge <- quantile(xys[, "y"], 0.95)
    list(EastEdge = xEdge, NorthEdge = yEdge)
  })
  edgesDataDT <- rbindlist(edgesDataCoords, idcol = "layerName")
  edgesDataDT[, type := "Data"]
  setnafill(edgesPredDT, type = "locf", cols = c("EastEdge", "NorthEdge"))

  edgesDT <- rbindlist(list(edgesDataDT, edgesPredDT))

  edgesDTLong <- melt(edgesDT, measure = patterns("Edge"))
  edgesDTLong[, value := c(NA, diff(value))/1e3, by = c("type", "variable")]
  edgesDTLong[, xy := gsub("Edge", "", variable)]
  Plots(edgesDTLong, plot_CentroidShift,
        title = "Predicted vs. Observed edge displacment each year",
        ggsaveArgs = list(width = 8, height = 10, units = "in", dpi = 300),
        filename = paste0(Par$type, ": Edge Displacement Pred vs Observed ", startToEnd))



  # Centroids
  centroidPredicted <- centroidChange(sim$predictedStack, sim$propPineRas)
  centroidPredicted[, variable := paste0("X", as.numeric(gsub("X", "", variable)) - 1)]
  setnames(centroidPredicted, old = c("yAtMax", "xAtMax"), new = c("yPredicted", "xPredicted"))

  centroidData <- centroidChange(sim$massAttacksStack, sim$propPineRas)
  setnames(centroidData, old = c("yAtMax", "xAtMax"), new = c("yData", "xData"))
  centroids <- centroidData[centroidPredicted]

  setnames(centroids, "variable", "layerName")
  centroidsLong <- melt(centroids, measure = patterns("Data|Pred"))
  centroidsLong[, type := gsub("x|y", "", variable)]
  centroidsLong[, xy := gsub("^(x|y).+", "\\1", variable)]
  centroidsLong[, xy := ifelse(xy == "x", "East", "North")]
  centroidsLong[, UTM := value]
  centroidsLong[, value := c(NA, diff(value))/1e3, by = c("variable")]

  Plots(centroidsLong, plot_CentroidShift, title = "Predicted vs. Observed centroid displacment each year",
        ggsaveArgs = list(width = 8, height = 10, units = "in", dpi = 300),
        filename = paste0(Par$type, ": Centroid Displacement Pred vs Observed ",
                                                     startToEnd))


  # Displacement Table
  centroids[, EastObs := c(NA, diff(xData)/1e3)]
  centroids[, NorthObs := c(NA, diff(yData)/1e3)]
  centroids[, EastPred := c(NA, diff(xPredicted)/1e3)]
  centroids[, NorthPred := c(NA, diff(yPredicted)/1e3)]
  centroids <- na.omit(centroids, cols = "EastPred")
  corEastDisplacement <- cor(centroids$EastObs, centroids$EastPred, method = "spearman")
  corNorthDisplacement <- cor(centroids$NorthObs, centroids$NorthPred, method = "spearman")
  sim$correlationsDisplacement <- c(East = corEastDisplacement, North = corNorthDisplacement)

  if (end(sim) - start(sim) >= 10) {
    centroids[, FiveYrPd := rep(c(paste0(layerName[1],"/",layerName[5]), paste0(layerName[6],"/",layerName[10])), each = 5)]

    centroids5Y <- centroids[, lapply(.SD, function(x) sum(x)), by = FiveYrPd,
                             .SDcols = c("EastObs", "NorthObs", "EastPred", "NorthPred")]
    centroids <- rbindlist(list(centroids, centroids5Y), fill = TRUE)
    centroids[is.na(layerName), layerName := FiveYrPd]
  }

  centroids2 <- as.data.table(t(centroids[, list(NorthObs, NorthPred, EastObs, EastPred)]))
  setnames(centroids2, centroids$layerName)
  centroids2 <- cbind(layerName = c("NorthObs", "NorthPred", "EastObs", "EastPred"), centroids2)


  sim$displacementTable <- centroids2 %>%
    gt() %>%
    tab_spanner(
      label = "Annual displacement (km/y)",
      columns = grep(value = TRUE, "^X.{4,4}$", colnames(centroids2))
    ) %>%
    tab_spanner(
      label = "Net displacement\n(km/5 y)",
      columns = grep(value = TRUE, "^X.+/.+", colnames(centroids2))
    )
  gtsave(sim$displacementTable,
         filename = file.path(outputPath(sim), "figures",
                              paste0(Par$type, ": displacementTable ", startToEnd, ".png")))


  return(invisible(sim))
}

##########################
plot_StudyAreaFn <- function(absk, cols, studyArea, mam, pred, propPineMap, minThreshold) {
  mam[mam[] == 0] <- NA
  pred[pred[] < log(minThreshold * 12)] <- NA
  pred[] <- pred[] * 0.7

  absksp <- sf::as_Spatial(absk)
  propPineMap[propPineMap[] <= 0.05] <- NA

  viewMode <- suppressMessages(identical(tmap_mode(), "view"))
  tmap_style("white")
  t1 <- tm_shape(absk) + tm_polygons(alpha = 0) + tm_graticules(alpha = 0.2) #tm_grid(projection = st_crs(4326), alpha = 0.15)
  t1 <- t1 + tm_shape(studyArea) + tm_polygons(alpha = 0, border.col = "lightblue")
  if (viewMode)
    t1 <- tm_basemap("OpenStreetMap") + t1 # only with tmap_mode("view")
  #t5 <- tm_tiles("OpenStreetMap")
  brks <- c(0, 6, 8, 10, 12, 14)
  t2 <- tm_shape(mam) + tm_raster(title = "Accumulated Damage", alpha = 0.8, palette = "YlOrRd",
                                  style = "fixed", breaks = brks, legend.show = FALSE)
  t3 <- tm_shape(propPineMap) + tm_raster(title = "Pine Cover (prop)", alpha = 0.6, palette = "Greens")
  t3b <- t3 + tm_legend(show = FALSE)
  t4 <- tm_shape(pred) + tm_raster(title = "MPB log(Attacked trees)", alpha = 0.8, palette = "YlOrRd", style = "fixed", breaks = brks)
  tpred <- t3 + t4 + t1 + tm_legend(show = TRUE, position = c("right", "top"))
  if (viewMode) {
    tpred
  } else {
    tpred <- tpred + tm_layout(title = "Predicted")
    tmam <- t3b + t2 + t1 + tm_compass(type  = "4star", position = c("right", "top"), north = 10) + tm_layout(title = "Observed")
    tmap_arrange(tmam, tpred)
  }
}



### plotting
plotFn <- function(sim) {
  pBest <- sim$fit_mpbSpreadOptimizer$optim$bestmem

  plotKernels(pBest)



  ## Histograms showing that there is no relationship between where beetle is and
  #  how much pine there is
  par(mfrow = c(3,4));
  for(i in 1:10/10) {
    cc <- sim$massAttacksStack$X2011[][sim$propPineRas[] < i & sim$propPineRas[] > (i - 0.1)];
    hist(log(cc), xlim = c(-1, 10), breaks = 15,
         main = paste0("Prop pine: ", i - 0.1, " to ", i), xlab = "log number attacked trees")
  }

  return(invisible(sim))
}

### spread
dispersalFit <- function(quotedSpread, propPineRas, studyArea, massAttacksDT, massAttacksStack,
                       rasterToMatch, maxDistance, # massAttacksRas,
                       params, #currentTime, bgSettlingProp,
                       type, reqdPkgs,
                       windDirStack, windSpeedStack, dispersalKernel) {

  # Make sure propPineRas is indeed a proportion
  mv <- maxValue(propPineRas)
  if (mv > 1)
    propPineRas[] <- propPineRas[] / 100 ## much faster than calc; ## TODO: allow different params per species

  # Set NAs on propPineRas to 0
  nas <- is.na(propPineRas[])
  propPineRas[nas] <- 0

  # Put objects in Global for objFun -- this is only when not using multi-machine cluster
  omitPastPines <- FALSE # was TRUE ... but maybe bad because there are many that have multiple years
  kern <- substr(tolower(dispersalKernel), start = 1, stop = 6)
  ll <- list()
  pars <- switch(kern,
         weibul = {
           ll$p <- do.call(c, params[c("p_meanDist", "p_advectionMag", "p_sdDist")])
           ll$lower <- c(5000,  0.01,  0.01)
           ll$upper <- c(105000, 2000, 3.5)
           ll
         },
         genera = {
           ll$p <- do.call(c, params[c("p_meanDist", "p_advectionMag", "p_sdDist", "p_nu")])#, "p_advectionDir")])
           ll$lower <- c(5000,  0.01,  0.3, 0.5)
           ll$upper <- c(105000, 2000, 1, 1)
           ll
         },
         expone = {
           ll$p <- do.call(c, params[c("p_meanDist", "p_advectionMag")])
           ll$lower <- c(5000,  0.1)
           ll$upper <- c(105000, 2000)
           ll
         },
         twodt = {
           ll$p <- do.call(c, params[c("p_meanDist", "p_advectionMag", "p_sdDist")])
           ll$lower <- c(5000,  0.1,  0.01)
           ll$upper <- c(105000, 2000, 0.5)
           ll
         },
  )
  p <- pars$p
  lower <- pars$lower
  upper <- pars$upper

  p_sdDist <- p["p_sdDist"] # sqrt(variance(of the Weibull distribution)) --> contributes to shape and scale parameters
  p_advectionMag <- p["p_advectionMag"]
  p[] <- sapply(seq_along(p), function(x) runif(1, lower[x], upper[x]))

  libPaths <- .libPaths()

  # Identify the objsToExport -- needs to be ALL the formalArgs of the objFun
  #   however, don't pass ones that are passed manually at the DEoptim call
  #   Generally, only pass "very small" objects in the DEoptim call. Large ones do this way (i.e.,
  #   export them first.
  objsToExport <- setdiff(formalArgs("objFun"),
                          # The following are small and passed at the time of the DEoptim call
                          c("p", "quotedSpread", "distanceFunction"))
  # need both the objsToExport vector, and the object called "objsToExport"
  objsToExport <- c(objsToExport, "reqdPkgs", "objsToExport", "libPaths")

  # put objsToExport into .GlobalEnv -- the DEoptim function gets all arguments from .GlobalEnv
  list2env(mget(objsToExport), envir = .GlobalEnv)

  # 213 appears slower than rest
  ipsNum <- c(189, 184, 217, 97, 106, 220, 213)
  ips <- paste0("10.20.0.", ipsNum)
  ips <- rep(ips, c(33, 16, 13, 23, 20, 5, 0))
  adjustments = 1#c(0.8, 1, 1, 1, 1, 1, 1)
  numCoresNeeded <- 110 # This is one per true core on 4 machines, 56, 56, 28, 28 cores each

  reqdPkgs <- grep(paste(collapse = "|", c("SpaDES.tools", "raster", "CircStats", "data.table",
                                           "purrr", "mpbutils", "gamlss")),
                   reqdPkgs, value = TRUE)
  stPre <- Sys.time()
  if (isTRUE(type %in% c("DEoptim", "fit"))) {
    cl <- clusterSetup(workers = ips, objsToExport = objsToExport,
                       reqdPkgs = reqdPkgs, libPaths = libPaths, doSpeedTest = 0, # fn = fn,
                       quotedExtra = quote(install.packages(c("rgdal", "rgeos", "sf",
                                                              "sp", "raster", "terra", "lwgeom"),
                                                             repos = "https://cran.rstudio.com")),
                       numCoresNeeded = numCoresNeeded, adjustments = adjustments)
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    message("Starting DEoptim")
    fit_mpbSpreadOptimizer <- DEoptim(fn = objFun,
                     lower = lower,#c(500, 1, -90, 0.9, 1.1, 5),
                     upper = upper,#c(30000, 10, 240, 1.8, 1.6, 30),
                     # reps = 1,
                     quotedSpread = quotedSpread,
                     control = DEoptim.control(cluster = cl,
                                               strategy = 2,
                                               itermax = 60,
                                               NP = numCoresNeeded))

    fit_mpbSpreadOptimizer <- colnamesToDEout(fit_mpbSpreadOptimizer, names(p))
    saveRDS(fit_mpbSpreadOptimizer, file = file.path("outputs", paste0("DEout_", format(stPre), ".rds")))
    sim2 <- get("sim", whereInStack("sim"))
    parms <- params(sim2)
    saveRDS(parms, file = file.path("outputs", paste0("parms_", format(stPre), ".rds")))
    try(parallel::stopCluster(cl), silent = TRUE)

  } else if (isTRUE(type == "runOnce")) {
    DEoutFileList <- dir("outputs", pattern = "fit_mpbSpreadOptimizer", full.names = TRUE)
    fit_mpbSpreadOptimizer <- readRDS(tail(DEoutFileList, 1))
    if (length(colnames(fit_mpbSpreadOptimizer$member$pop)) == 0) {
      message("Please add names to parameter vector in fit_mpbSpreadOptimizer")
      fit_mpbSpreadOptimizer <- colnamesToDEout(fit_mpbSpreadOptimizer, c("p_meanDist", "p_advectionMag", "p_sdDist", "p_rescaler"))
    }
    pBest <- fit_mpbSpreadOptimizer$optim$bestmem[] # use [] to keep names
    pBest["p_sdDist"] <- 0.08
    pBest <- p
    out <- objFun(quotedSpread = quotedSpread, # reps = 1,
                  p = pBest,
                  massAttacksStack = massAttacksStack, massAttacksDT = massAttacksDT)
  } else if (isTRUE(type == "optim")) {
    # stop("This has not been maintained and appears to be not capable of estimating parameters")
    numCoresNeeded <- 5
    cl <- clusterSetup(rep("localhost", numCoresNeeded), objsToExport = objsToExport,
                       reqdPkgs = reqdPkgs, libPaths = libPaths, doSpeedTest = FALSE
                        )
    on.exit(try(parallel::stopCluster(cl)), add = TRUE)
    fit_mpbSpreadOptimizer <- optimParallel::optimParallel(fn = objFun, par = p,
                                             lower = lower,
                                             upper = upper,
                                             method = "L-BFGS-B",
                                             quotedSpread = quotedSpread,
                                             control = list(trace = 3, factr = 1e6),
                                             parallel = list(cl = cl, forward = TRUE, loginfo  = TRUE)
                                             )
    saveRDS(fit_mpbSpreadOptimizer, file = file.path("outputs", paste0("optimOut_", format(stPre), ".rds")))
  }
  try(parallel::stopCluster(cl))
  stPost <- Sys.time()
  print(difftime(stPost, stPre))

  return(invisible(fit_mpbSpreadOptimizer))
}



spiralDistances <- function(pixelGroupMap, maxDis, cellSize) {
  spiral <- which(focalWeight(pixelGroupMap, maxDis, type = "circle")>0, arr.ind = TRUE) -
    ceiling(maxDis/cellSize) - 1
  spiral <- cbind(spiral, dists = sqrt( (0 - spiral[,1]) ^ 2 + (0 - spiral[, 2]) ^ 2))
  spiral <- spiral[order(spiral[, "dists"], apply(abs(spiral), 1, sum),
                         abs(spiral[, 1]), abs(spiral[, 2])),, drop = FALSE]
}

distanceFunction <- function(dist, angle, landscape, fromCell, toCells, # nextYrVec,
                             propPineRas, p_advectionMag, windDir, p_meanDist, p_sdDist,
                             dispersalKernel, maxDistance, windSpeed, p_nu = NULL) {

  if (p_advectionMag > 0) {
    fromXY <- xyFromCell(propPineRas, fromCell)
    toXYs <-  xyFromCell(propPineRas, toCells)
    if (length(windDir) > 1) {
      windDir <- windDir[fromCell]
      windSpeed <- windSpeed[fromCell]
    }
    advectionXY <- c(x = sin(rad(windDir))*p_advectionMag*windSpeed, y = cos(rad(windDir))*p_advectionMag*windSpeed)
    x <- toXYs[, 'x'] - advectionXY[1]/p_meanDist*dist
    y <- toXYs[, 'y'] - advectionXY[2]/p_meanDist*dist
    newTo <- cbind(x, y)
    colnames(newTo) <- c("x", "y")
    if (!is.null(dim(dist))) dist <- dist[,1, drop = TRUE]
    dists <- cbind(origX = toXYs[, "x"], origY = toXYs[, "y"],
                   origDist = dist,
                   .pointDistance(from = fromXY, newTo, angles = FALSE))
    dist <- dists[, "dists"]
  }

  prob <- kernelFn(p_meanDist, p_sdDist, dispersalKernel, dist, p_nu = p_nu)

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


grow <- function(massAttacksDT, dataset, massAttacksStack, growthData, year = NULL) {
  ## determine the actual growth based on the actual number of attacked trees/ha
  layerNames <- unique(massAttacksDT$layerName)
  if (!is.null(year)) layerNames <- grep(year, layerNames, value = TRUE)
  massAttacksDT[layerName %in% layerNames, greenTreesYr_t := xt(ATKTREES, CLIMATE, dataset,
                                       massAttacksStack, growthData)]

  return(massAttacksDT)
}


hist.DEoptim <- function(DEobj, paramNames) {
  dt <- as.data.table(DEobj$member$pop)
  if (missing(paramNames))
    paramNames <- names(MPBfit$fit_mpbSpreadOptimizer$optim$bestmem)
  colnames(dt) <- paramNames
  dt <- melt(dt, measure.vars = seq(NCOL(dt)))
  gg <- ggplot(dt, aes(x = value)) +
    geom_histogram(bins = 15) +
    facet_grid(. ~ variable, scales = "free") +
    theme_bw()

  gg
}

summary2.DEoptim <- function(DEobj, title) {
  if (!missing(title))
    cat(paste(title, "\n"))
  cat(paste("Likelihood: ", round(DEobj$optim$bestval, 0), "\n"))
  cat(paste("Mean params: ", paste(round(apply(DEobj$member$pop, 2, mean), 2), collapse = ", "), "\n"))
  cat(paste("Best params: ", paste(round(DEobj$optim$bestmem, 2), collapse = ", "), "\n"))
}

plotKernels <- function(p, reps = 10) {
  if ( (is.na(p["meanDistSD"])))
    reps = 1
  nrows <- max(1, floor(sqrt(reps + 1)))
  ncols <- max(1, ceiling((reps + 1)/nrows))
  par(mfrow = c(nrows, ncols))
  if (is.null(names(p))) stop("p must be named and have these 3: p_meanDist, p_sdDist, meanDistSD")
  for (i in 1:reps) {
    p_meanDist <- p["p_meanDist"]
    if (!is.na(p["meanDistSD"]))
      p_meanDist <- rlnorm(1, log(p["p_meanDist"]), log(p["meanDistSD"]))

    p_sdDist = p["p_sdDist"]
    shapeScale <- weibullShapeScale(p_meanDist = p_meanDist, p_sdDist = p_sdDist)
    ww <- rweibull(1e4, shape = shapeScale$shape, scale = shapeScale$scale)
    hist(ww, xlim = c(0, 14e4), main = "", xlab = "Distance (m) from source")
    abline(v = p_meanDist, col = "red")
    abline(v = p[["p_meanDist"]], col = "green")
  }
  plot(type = "n", x = 1:10, axes = F, xlab = "", ylab = "")
  legend(x = 2, y = 8, col = c("red", "green"), lty= 1,
         legend = c("annual average", "global average"),
         cex = 1.5)
  mtext(text = "Dispersal kernels", outer = TRUE, line = -3, xpd = TRUE)
}

kernelFn <- function(p_meanDist, p_sdDist, dispersalKernel, dist, p_nu) {
  dispersalKernel <- tolower(dispersalKernel)
  if (grepl("weib", dispersalKernel)) {
    shapeScale <- weibullShapeScale(p_meanDist, p_sdDist)
    prob <- dweibull(dist, scale = shapeScale$scale, shape = shapeScale$shape)
  } else if (grepl("twodt", dispersalKernel)) {
    muConvert <- mpbutils:::kernel_twoDT_mean(1, 1)
    mu <- p_meanDist / muConvert
    prob <- kernel_twoDT(dist, mu, p = p_sdDist)
  } else if (grepl("generalgamma", dispersalKernel)) {
    prob <- dGG(dist, mu=p_meanDist, sigma=p_sdDist, nu=p_nu, log = F)
  } else {
    # Using exponential
    prob <- dexp(dist, rate = 1/p_meanDist)
  }
  nonFinites <- !is.finite(prob) # captures NaN and -Inf -- different ways that these d* function deal with zeros
  if (any(nonFinites))
    prob[nonFinites] <- 0
  prob
}

weibullShapeScale <- function(p_meanDist, p_sdDist) {
  mn <- (p_meanDist)
  sd <- mn/p_sdDist # 0.8 to 2.0 range
  shape <- (sd/mn)^(-1.086)
  scale <- mn/exp(lgamma(1+1/shape))
  return(list(scale = scale, shape = shape))
}


stackFromDT <- function(predictedDT, rasterTemplate, colForStack = "ATKTREES") {
  yrNames <- unique(predictedDT$layerName)
  stackDataAsMatrix <- matrix(rep(NA, length(yrNames) * ncell(rasterTemplate)),
                    ncol = length(yrNames))
  yrNamesPlus1 <- yrNamesPlus1(yrNames)
  colnames(stackDataAsMatrix) <- yrNamesPlus1

  cols <- c("layerName", "pixel", colForStack)
  outW <- dcast(predictedDT[, ..cols], pixel  ~ layerName, value.var = colForStack)
  cols <- grep("^X", colnames(outW), value = TRUE)
  gg <- rep(outW$pixel, NCOL(outW[, ..cols]))
  ff <- rep(as.integer(factor(colnames(stackDataAsMatrix))), each = NROW(outW))
  stackDataAsMatrix[cbind(gg, ff)] <- as.matrix(outW[, ..cols])

  rasterTemplate[] <- NA
  predictedStack <- raster::stack(lapply(yrNames, function(lay) rasterTemplate))
  predictedStack[] <- stackDataAsMatrix
  predictedStack <- raster::stack(predictedStack)
  predictedStack
}

yrNamesPlus1 <- function(yrNames) {
  yrNames <- gsub("X([[:digit:]]{4,4})", "\\1", yrNames)
  paste0("X", as.integer(yrNames) + 1)
}

yrNamesMinus1 <- function(yrNames) {
  yrNames <- gsub("X([[:digit:]]{4,4})", "\\1", yrNames)
  paste0("X", as.integer(yrNames) - 1)
}

createStackPredAndObs <- function(massAttacksStack, predictedStack, threshold = 1) {
  yrNamesPlus1 <- names(predictedStack)
  yrNames <- yrNamesMinus1(yrNamesPlus1)
  stk <- raster::stack()
  for (lay in seq(yrNamesPlus1)) {
    setColors(predictedStack[[yrNamesPlus1[lay]]], 8) <- "Reds"
    AllThree <- raster(predictedStack[[yrNamesPlus1[lay]]])
    AllThree[] <- 0
    if (yrNames[lay] %in% names(massAttacksStack)) {
      whCA <- which(!is.na(massAttacksStack[[lay]][]))
      AllThree[whCA] <- 1
    }
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
  raster::stack(stk)
}


stacksPredVObs <- function(massAttacksStack, predictedStack, propPineMap, threshold = 10, plot.it = FALSE) {
  yrNamesPlus1 <- names(predictedStack)
  yrNames <- yrNamesMinus1(yrNamesPlus1)
  # yrNames <- names(massAttacksStack)
  # yrNamesPlus1 <- yrNamesPlus1(yrNames)
  stk <- raster::stack()
  ROCs <- list()
  if (isTRUE(plot.it)) {
    numR <- floor(sqrt(length(yrNamesPlus1)))
    par(mfrow = c(numR, ceiling(length(yrNamesPlus1)/numR)))
  }
  for (lay in seq(yrNamesPlus1)) {

    preds <- predictedStack[[yrNamesPlus1[lay]]][]
    whHasVal <- if (yrNames[lay] %in% names(massAttacksStack)) {
      which(!is.na(preds) |
                          !is.na(massAttacksStack[[lay]][]) |
                          propPineMap[] > 0)
    } else {
      which(!is.na(preds))
    }
    preds[is.na(preds)] <- 0

    # setColors(predictedStack[[yrNamesPlus1[lay]]], 8) <- "Reds"
    Both <- raster(predictedStack[[yrNamesPlus1[lay]]])
    best <- raster(Both)
    # whCA <- which(!is.na(massAttacksStack[[yrNamesPlus1[lay]]][]))
    Both[] <- 0
    # Both[whCA] <- 1
    if (yrNamesPlus1[lay] %in% names(massAttacksStack)) {
      whNY <- which(!is.na(massAttacksStack[[yrNamesPlus1[lay]]][]))
      Both[whNY] <- Both[whNY] + 1
    }
    whPM <- which(!is.na(predictedStack[[yrNamesPlus1[lay]]][]) &
                    predictedStack[[yrNamesPlus1[lay]]][] > 0)
    Both[whPM] <- Both[whPM] + 2
    Both[propPineMap[] > 0 & (is.na(Both[]) | Both[] == 0) ] <- 4
    levels(Both) <- data.frame(ID = 1:4,
                               Label = c("Obs", "Pred", "Pred&Obs", "Neither"))
    setColors(Both, n = 4) <- "Set2"
    Both@legend@colortable[1:4] <- c("yellow", "blue", "green", "purple")

    dt <- data.table(pred = preds[whHasVal], outcome = Both[][whHasVal])
    dt[, Attacked := ifelse(outcome %in% c(1, 3), 1, 0)]
    dt <- setorderv(dt, c("pred", "Attacked"), order = c(-1, -1))
    dt <- dt[pred > 0]
    dt[, AttackedSum := cumsum(Attacked)]
    dt[, NoAttackSum := cumsum(Attacked == 0)]
    if (!all(dt$Attacked == 0)) {
      dt[, Sens := AttackedSum / max(AttackedSum, na.rm = TRUE)]
      dt[, OneMinusSpec := NoAttackSum / max(NoAttackSum)]
      ROCs[[lay]] <- data.table(year = yrNamesPlus1[lay],
                                bestAnnualThreshold = dt$pred[which.max(dt$Sens + (1 - dt$OneMinusSpec))],
                                threshold = threshold,
                                sumSensSpecAtThreshold = (dt$Sens + (1 - dt$OneMinusSpec))[which(dt$pred < threshold)[1]],
                                sumSensSpecAtBest = max(dt$Sens + (1 - dt$OneMinusSpec)),
                                ROC = sum(dt$Sens)/NROW(dt))
      if (isTRUE(plot.it)) {
        plot(dt$OneMinusSpec, dt$Sens, type = "l", main = yrNamesPlus1[lay])
        abline(a = 0, b = 1, lty = "dashed")
        # print(paste("year: ", yrNamesPlus1[lay],", threshold: ", threshold, ", ROC: ", round(ROCs[[lay]]$ROC, 2)))
      }
    }
    stk <- addLayer(stk, Both)

  }
  names(stk) <- paste0("pred", yrNamesPlus1)
  stk[which(stk[] == 0)] <- NA
  list(stk = raster::stack(stk), ROCs = rbindlist(ROCs))
}


centroidChange <- function(stk, propPineRas) {
  mats <- xyFromCell(stk, cell = seq(ncell(stk)))
  masDT <- as.data.table(stk[])
  bb <- data.table(propPine = propPineRas[], mats, masDT, pixel = seq(ncell(stk)))
  cc <- melt(bb, measure.vars = patterns("^X"), value.name = "NumAttackedTrees")
  cc <- cc[!is.na(NumAttackedTrees)]
  cc[, `:=`(sumByX = sum(NumAttackedTrees, na.rm = TRUE)),
     by = c("x", "variable")]
  cc[, `:=`(sumByY = sum(NumAttackedTrees, na.rm = TRUE)),
     by = c("y", "variable")]
  xs <- cc[, .SD[1], by = c("x", "variable")]
  ys <- cc[, .SD[1], by = c("y", "variable")]
  setkeyv(xs, c("variable", "x", "y"))#; setkey(uniX, "x")
  setkeyv(ys, c("variable", "x", "y"))#; setkey(uniX, "x")
  keepsY <- c("variable", "propPine", "x", "sumByX")
  keepsX <- c("variable", "propPine", "y", "sumByY")
  setkeyv(xs, c("variable", "y"))#; setkey(uniX, "x")
  setkeyv(ys, c("variable", "x"))#; setkey(uniX, "x")
  ys <- unique(ys[, ..keepsY], by = c("x", "variable"))
  xs <- unique(xs[, ..keepsX], by = c("y", "variable"))
  ys[, cumsumSumByX := cumsum(sumByX), by = c("variable")]
  xs[, cumsumSumByY := cumsum(sumByY), by = c("variable")]
  yAtMax <- xs[, list(yAtMax = y[max(which(cumsumSumByY < (max(cumsumSumByY)/2) ))]), by = "variable"]
  xAtMax <- ys[, list(xAtMax = x[max(which(cumsumSumByX < (max(cumsumSumByX)/2) ))]), by = "variable"]
  centroids <- yAtMax[xAtMax]
  centroids
}

plot_CentroidShift <- function(centroids, title) {
  ggplot(centroids[], aes(x = layerName, y = value, group = xy, color = type)) +
    geom_point() +
    facet_grid(vars(xy), scales = "free_y") +
    ggplot2::theme_bw() +
    ggtitle(title)
}

plot_HistsOfDEoptimPars <- function(fit_mpbSpreadOptimizer, title) {
  (paramHists <- fit_mpbSpreadOptimizer %>%
    ggplot( aes(x=value, fill=Parameter)) +
    geom_histogram() +# color="#e9ecef", alpha=0.6, position = 'identity') +
    facet_grid(. ~ Parameter, scales = "free") +
    theme_bw() +
    ggtitle(title))
}

growPredict <- function(sim) {
  if (time(sim) == start(sim)) {
    massAttacksDTforPredict <- data.table::copy(sim$massAttacksDT)
    massAttacksDTforPredict <- massAttacksDTforPredict[grep(time(sim), layerName), ]
  } else {
    massAttacksDTforPredict <- data.table::copy(sim$predictedDT)
  }
  if ("predYear" %in% colnames(massAttacksDTforPredict)) {
    massAttacksDTforPredict[, layerName := predYear]
  }

  massAttacksDTforPredict <- massAttacksDTforPredict[grep(time(sim), layerName)] # remove all future data
  grow(massAttacksDTforPredict, P(sim)$dataset,
         sim$massAttacksStack, mod$growthData, year = time(sim))
}

colnamesToDEout <- function(fit_mpbSpreadOptimizer, paramNames = c("p_meanDist", "p_advectionMag", "p_sdDist", "p_rescaler")) {

  if (length(colnames(fit_mpbSpreadOptimizer$member$pop)) == 0) {
    message(crayon::green("Please ensure parameter vector in fit_mpbSpreadOptimizer is named: ", paste(collapse = ", ", paramNames)))
    names(fit_mpbSpreadOptimizer$optim$bestmem) <- paramNames
    names(fit_mpbSpreadOptimizer$member$lower) <- paramNames
    names(fit_mpbSpreadOptimizer$member$upper) <- paramNames
    colnames(fit_mpbSpreadOptimizer$member$bestmemit) <- paramNames
    colnames(fit_mpbSpreadOptimizer$member$pop) <- paramNames
  }

  fit_mpbSpreadOptimizer
}

estimateMinAttackThresh <- function(p_minDensityThresh, meanAttackStk, predictedStk, propPineRast) {
  stackPredVObs <- stacksPredVObs(meanAttackStk, predictedStk, propPineRast, threshold = p_minDensityThresh[[1]])
  2 - mean(stackPredVObs$ROCs$sumSensSpecAtThreshold)
}


makeAnnualMeansDT <- function(pred, obs) {
  predM <- pred[]
  predMean <- colMeans(predM, na.rm = TRUE)
  obsM <- obs[]
  obsMean <- colMeans(obsM, na.rm = TRUE)
  dt <- data.table(`Predicted Mean` = predMean,
                   `Observed Mean` = obsMean, Year = names(predMean))
  dt
}
plot_ObsVPredAbund1 <- function(dt) {
  ggplot(dt, aes(x = `Observed Mean`, y = `Predicted Mean`)) +
    geom_point() +
    theme_bw() +
    xlab("Observed mean attack density, per year") +
    ylab("Predicted mean attack density, per year")
}

plot_ObsVPredAbund2 <- function(dt) {
  dtL <- melt(dt, id.vars = "Year")
  ggplot(dtL, aes(x = Year, y = value, group = variable, col = variable)) +
    geom_line() +
    theme_bw() +
    scale_y_continuous(trans='log') +
    xlab("Year") +
    ylab("Predicted and Observed mean attack density, per year")
}
