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
                  "CircStats", "data.table", "DEoptim", "EnvStats", "ggplot2",
                  "PredictiveEcology/SpaDES.core@modifyList2 (>= 1.0.8.9004)",
                  "PredictiveEcology/LandR@cluster (>= 1.0.4.9019)",
                  "parallelly",
                  "PredictiveEcology/pemisc@development (>= 0.0.3.9001)",
                  "PredictiveEcology/mpbutils (>= 0.1.2)", "purrr",
                  "quickPlot", "raster", "RColorBrewer",
                  "PredictiveEcology/reproducible@DotsBugFix (>= 1.2.7.9009)",
                  "PredictiveEcology/SpaDES.tools@development (>= 0.3.7.9021)"),
  parameters = rbind(
    defineParameter("dataset", "character", "Boone2011", NA, NA, "Which dataset to use for stand dynamic model fitting. One of 'Boone2011' (default), 'Berryman1979_fit', or 'Berryman1979_forced'. Others to be implemented later."),
    defineParameter("growthInterval", "numeric", 1, NA, NA, "This describes the interval time between growth events"),
    defineParameter("advectionDir", "numeric", 90, 0, 359.9999,
                    "The direction of the spread bias, in degrees from north"),
    defineParameter("advectionMag", "numeric", 5e3, NA, NA,
                    "The magnitude of the wind effect on spread. This number is multiplied by the wind speed (which averages 9.6 km/h in the historical dataset)"),
    defineParameter("bgSettlingProp", "numeric", 0.1, 0, 1,
                    "The proportion of beetles that settle from those that could potentially settle, even if no pine"),
    defineParameter("maxDistance", "numeric", 1.4e5, NA, NA,
                    "The maximum distance to allow for pair-wise from-to dispersal pairs"),
    defineParameter("meanDist", "numeric", 1e4, NA, NA,
                    "Expected dispersal distance (m); ~63% go less than this distance"),
    defineParameter("quotedSpread", "language", NULL, NA, NA,
                    "A spread function that is contained within a quote(). If left at NULL or NA, it will use the module version"),
    defineParameter("sdDist", "numeric", 1.2, NA, NA,
                    paste0("The dispersion term for the Weibull dispersal kernel: contributes to shape and scale parameters;",
                    "sqrt(variance(of the Weibull distribution)) ")),
    defineParameter("type", "character", "DEoptim", NA, NA,
                    "One of three modes of running this module: DEoptim, optim, runOnce or nofit"),
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
    createsOutput("pBest", "numeric", "Best p values coming from DEoptim"),
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
           sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "dispersal", eventPriority = 5.5)
           sim <- scheduleEvent(sim, time(sim), "mpbRedTopSpread", "predict", 6.5)
           # sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "explore", 7.5)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "mpbRedTopSpread", "plot", 7.5)
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
           # The outputs for this are saved to disk ... so no "return" sim
           sim$pBest <- dispersal2(quotedSpread = P(sim)$quotedSpread,
                             propPineRas = sim$propPineRas, studyArea = sim$studyArea,
                             massAttacksDT = sim$massAttacksDT,
                             massAttacksStack = sim$massAttacksStack,
                             maxDistance = P(sim)$maxDistance,
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

           # sim <- dispersal(sim)
           # sim <- scheduleEvent(sim, time(sim) + 1, "mpbRedTopSpread", "dispersal", eventPriority = 5.5)
         },
         "predict" = {
           Predict(sim)
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

  if (is.null(P(sim)$quotedSpread))
    params(sim)[[currentModule(sim)]]$quotedSpread <- quotedSpreadDefault


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

Predict <- function(sim) {
  objs <- dir("outputs", pattern = "DEout", full.names = TRUE)
  DEout <- readRDS(tail(objs, 1))
  DEout <- readRDS(tail(objs, 2)[1])

  p <- sim$pBest # gets the parameter names
  p[] <- DEout$optim$bestmem # replace with values from disk object


  DEout <- as.data.table(DEout$member$pop)
  setnames(DEout, new = names(p))
  DEout <- melt(DEout, measure = colnames(DEout), variable.name = "Parameter")

  Plots(DEout, plotHistsOfDEoptimPars, title = "Parameter histograms from DEoptim",
        filename = paste0("Histograms of parameters from DEoptim ", Sys.time()))


  # Way to identify quantiles
  if (FALSE) {
    ll <- weibullShapeScale(p['meanDist'], p['sdDist']);
    qweibull(0.99, shape = ll$shape, scale = ll$scale)
  }

  #########################################################
  #################### PREDICT MAPS
  #########################################################

  opts2 <- options(reproducible.cacheSpeed = "fast")
  on.exit(options(opts2))
  # Browse[1]> fastdigest(quotedSpread)
  # [1] "6ee73bef9a187c2404d856d22b7a426a"
  sim$predictedDT <- Cache(predictQuotedSpread,
                       massAttacksDT = sim$massAttacksDT,
                       massAttacksStack = sim$massAttacksStack,
                       windDirStack = sim$windDirStack,
                       windSpeedStack = sim$windSpeedStack,
                       propPineRas = sim$propPineRas,
                       pineThreshold = 0.3,
                       dispersalKernel = "Weibull",
                       maxDistance = P(sim)$maxDistance,
                       quotedSpread = P(sim)$quotedSpread, # doesn't cache correctly
                       .cacheExtra = format(P(sim)$quotedSpread), # cache this instead
                       p = p,
                       omitArgs = "quotedSpread",
                       subsampleFrom = subsampleFrom
  )
  sim$predictedStack <- Cache(stackFromDT, sim$massAttacksStack, sim$predictedDT)
  #bb <- predictedStack$X2011[][propPineRas[] < 0.3]
  #cc <- massAttacksStack$X2011[][propPineRas[] < 0.1]

  lapply(seq(nlayers(sim$massAttacksStack)), function(x)
    cor.test(sim$predictedStack[][, x], sim$massAttacksStack[][, x], use = "complete.obs"))

  corByYr <- cor(sim$predictedStack[], sim$massAttacksStack[], use = "complete.obs")
  (corByYrAvg <- mean(diag(corByYr)))

  plotStack <- Cache(stacksForPlot, sim$massAttacksStack, sim$predictedStack, threshold = 55)

  fnHash <- fastdigest(format(stacksPredVObs))
  fn <- function(p, mas, ps, pp) {
    stackPredVObs <- stacksPredVObs(mas, ps, pp, threshold = p[[1]])
    2 - mean(stackPredVObs$ROCs$sumSensSpecAtThreshold)
    # 1 - mean(stackPredVObs$ROCs$ROC)
  }
  p <- 10

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
  ss <- sum(mam)
  ss[] <- log(ss[] + 1)
  # 3setColors(ss, n = 5) <- c("yellow", "orange" , "red", "purple", "blue")

  # cumulative maps
  pred <- sim$predictedStack
  for (lay in names(pred)) {
    pred[[lay]][is.na(pred[[lay]][])] <- 0
  }
  pp <- sum(pred)
  pp[] <- log(pp[] + 1)


  ##########################
  ggplotStudyAreaFn <- function(absk, cols, studyArea, mam, pred, propPineMap) {
    mam[mam[] == 0] <- NA
    pred[pred[] < log(60 * 12)] <- NA
    pred[] <- pred[] * 0.7

    absksp <- sf::as_Spatial(absk)
    propPineMap[propPineMap[] <= 0.05] <- NA
    browser()

    current.mode <- tmap_mode("plot")
    # current.mode <- tmap_mode("view")

    #t0 <- tm_basemap("OpenStreetMap") # only with tmap_mode("view")
    #t5 <- tm_tiles("OpenStreetMap")
    tmap_style("white")
    brks <- c(0, 6, 8, 10, 12, 14)
    t1 <- tm_shape(absk) + tm_polygons(alpha = 0)
    t2 <- tm_shape(mam) + tm_raster(title = "Accumulated Damage", alpha = 0.8, palette = "YlOrRd",
                                    style = "fixed", breaks = brks, legend.show = FALSE)
    t3 <- tm_shape(propPineMap) + tm_raster(title = "Pine Cover (prop)", alpha = 0.6, palette = "Greens")
    t3b <- t3 + tm_legend(show = FALSE)
    t4 <- tm_shape(pred) + tm_raster(title = "MPB log(Attacked trees)", alpha = 0.8, palette = "YlOrRd", style = "fixed", breaks = brks)
    tmam <- t3b + t2 + t1 + tm_layout(title = "Observed")
    tpred <- t3 + t4 + t1  + tm_compass(type  = "4star", position = c("right", "top"), north = 10) + tm_layout(title = "Predicted")
    tmap_arrange(tmam, tpred)
    tmam + tpred
  }

  browser()
  ggplotStudyAreaFn(sim$absk, "darkgreen", sf::st_as_sf(sim$studyArea), ss,
                    pred = pp,
                    propPineMap = sim$propPineRas)

  system.time(optimOut <- Cache(optimize, f = fn, interval = c(0.001, 200),
                                mas = sim$massAttacksStack, .cacheExtra = fnHash,
                                ps = sim$predictedStack, pp = sim$propPineRas,
                                tol = 0.001
  ))

  thresholdBest <- optimOut$minimum
  browser()
  stackPredVObs <- stacksPredVObs(sim$massAttacksStack, sim$predictedStack,
                                  apmBuff, threshold = thresholdBest,
                                  plot = TRUE)



  mean(stackPredVObs$ROCs$ROC)
  stackPredVObs$ROCs

  X2011 <- sim$predictedStack$X2011
  X2011gtThresh <- X2011
  X2011pred <- X2011
  X2011gtThresh[X2011 < thresholdBest] <- NA
  X2011pred[sim$massAttacksStack$X2011[] > 0]
  p0 <- levelplot(X2011, par.settings = GrTheme(brewer.pal(9, "Greys")))
  p1 <- levelplot(X2011gtThresh, par.settings = RdBuTheme(region = brewer.pal(9, 'Greens')))
  p2 <- levelplot(X2011gtThresh, par.settings = RdBuTheme(region = brewer.pal(9, 'Greens')))
  p0 + p1

  gg <- gplot(sim$predictedStack$X2011)




  Plot(stackPredVObs$stk)

  browser()

  ss <- lapply(raster::unstack(stackPredVObs), function(ras) {
    freqs <- table(ras[])
    sensitivity <- freqs[3]/(sum(freqs[3], freqs[1]))
    specificity <- freqs[4]/(sum(freqs[4], freqs[2]))
    return(list(sensitivity = sensitivity, specificity = specificity))
  })


  browser()

  # Edges
  edgesPredicted <- sim$predictedStack
  edgesPredicted[edgesPredicted < thresholdBest] <- NA
  edgesPredicted <- raster::stack(edgesPredicted)
  edgesList <- raster::unstack(edgesPredicted)
  names(edgesList) <- names(edgesPredicted)
  edgesCoords <- lapply(edgesList, function(edge) {
    xys <- xyFromCell(edge, which(!is.na(edge[])))
    xEdge <- quantile(xys[, "x"], 0.95)
    yEdge <- quantile(xys[, "y"], 0.95)
    list(EastEdge = xEdge, NorthEdge = yEdge)
  })
  edgesPredDT <- rbindlist(edgesCoords, idcol = "Year")
  edgesPredDT[, type := "Prediction"]

  edgesData <- sim$massAttacksStack
  edgesData <- raster::stack(edgesData)
  edgesDataList <- raster::unstack(edgesData)
  names(edgesDataList) <- names(edgesData)
  edgesDataCoords <- lapply(edgesDataList, function(edge) {
    xys <- xyFromCell(edge, which(!is.na(edge[])))
    xEdge <- quantile(xys[, "x"], 0.95)
    yEdge <- quantile(xys[, "y"], 0.95)
    list(EastEdge = xEdge, NorthEdge = yEdge)
  })
  edgesDataDT <- rbindlist(edgesDataCoords, idcol = "Year")
  edgesDataDT[, type := "Data"]
  setnafill(edgesPredDT, type = "locf", cols = c("EastEdge", "NorthEdge"))

  edgesDT <- rbindlist(list(edgesDataDT, edgesPredDT))

  edgesDTLong <- melt(edgesDT, measure = patterns("Edge"))
  edgesDTLong[, value := c(NA, diff(value))/1e3, by = c("variable")]
  edgesDTLong[, xy := gsub("Edge", "", variable)]
  Plots(edgesDTLong, plotCentroidShift, title = "Predicted vs. Observed edge displacment each year")



  # Centroids
  centroidPredicted <- centroidChange(sim$predictedStack, propPineRas)
  centroidPredicted[, variable := paste0("X", as.numeric(gsub("X", "", variable)) - 1)]
  setnames(centroidPredicted, old = c("yAtMax", "xAtMax"), new = c("yPredicted", "xPredicted"))

  centroidData <- centroidChange(sim$massAttacksStack, propPineRas)
  setnames(centroidData, old = c("yAtMax", "xAtMax"), new = c("yData", "xData"))
  centroids <- centroidData[centroidPredicted]

  setnames(centroids, "variable", "Year")
  centroidsLong <- melt(centroids, measure = patterns("Data|Pred"))
  centroidsLong[, type := gsub("x|y", "", variable)]
  centroidsLong[, xy := gsub("^(x|y).+", "\\1", variable)]
  centroidsLong[, xy := ifelse(xy == "x", "East", "North")]
  centroidsLong[, UTM := value]
  centroidsLong[, value := c(NA, diff(value))/1e3, by = c("variable")]

  Plots(centroidsLong, plotCentroidShift, title = "Predicted vs. Observed centroid displacment each year")

  centroids[, EastObs := c(NA, diff(xData)/1e3)]
  centroids[, NorthObs := c(NA, diff(yData)/1e3)]
  centroids[, EastPred := c(NA, diff(xPredicted)/1e3)]
  centroids[, NorthPred := c(NA, diff(yPredicted)/1e3)]
  centroids <- na.omit(centroids, cols = "EastPred")
  centroids[, FiveYrPd := rep(c(paste0(Year[1],"/",Year[5]), paste0(Year[6],"/",Year[10])), each = 5)]

  centroids5Y <- centroids[, lapply(.SD, function(x) sum(x)), by = FiveYrPd,
                           .SDcols = c("EastObs", "NorthObs", "EastPred", "NorthPred")]
  centroids <- rbindlist(list(centroids, centroids5Y), fill = TRUE)
  centroids[is.na(Year), Year := FiveYrPd]

  centroids2 <- as.data.table(t(centroids[, list(NorthObs, NorthPred, EastObs, EastPred)]))
  setnames(centroids2, centroids$Year)
  centroids2 <- cbind(Year = c("NorthObs", "NorthPred", "EastObs", "EastPred"), centroids2)
  displacementTable <- centroids2 %>%
    gt() %>%
    tab_spanner(
      label = "Annual displacement (km/y)",
      columns = grep(value = TRUE, "^X.{4,4}$", colnames(centroids2))
    ) %>%
    tab_spanner(
      label = "Net displacement\n(km/5 y)",
      columns = grep(value = TRUE, "^X.+/.+", colnames(centroids2))
    )
  gtsave(displacementTable, filename = file.path(outputPath(sim), "displacementTable.png"))


  return(sim)


}
### plotting
plotFn <- function(sim) {
  p <- sim$pBest

  plotKernels(p)


  par(mfrow = c(1,2))
  plot(colMeans(sim$climateSuitabilityMaps[[names(sim$massAttacksStack)]][], na.rm = TRUE))
  plot(colMeans(sim$massAttacksStack[[names(sim$massAttacksStack)]][], na.rm = TRUE))

  ## Histograms showing that there is no relationship between where beetle is and
  #  how much pine there is
  par(mfrow = c(4,6));
  for(i in 1:10/10) {
    cc <- massAttacksStack$X2011[][propPineRas[] < i & propPineRas[] > (i - 0.1)];
    hist(log(cc), xlim = c(-1, 10), breaks = 15,
         main = paste0("Prop pine: ", i - 0.1, " to ", i), xlab = "log number attacked trees")
  }


  r <- raster(sim$massAttacksStack[[1]]) ## use a template
  # atkMap <- amc::dt2raster(sim$massAttacksDT, r, "ATKTREES")
  # setColors(atkMap) <- RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")
  # Plot(atkMap, title = "Simulated Attacks")

  return(invisible(sim))
}

### spread
dispersal2 <- function(quotedSpread, propPineRas, studyArea, massAttacksDT, massAttacksStack,
                       rasterToMatch, maxDistance, # massAttacksRas,
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
  mv <- maxValue(propPineRas)
  #
  if (mv > 1)
    propPineRas[] <- propPineRas[] / 100 ## much faster than calc; ## TODO: allow different params per species

  nas <- is.na(propPineRas[]) # & !is.na(rasterToMatch[])
  propPineRas[nas] <- 0

  # Put objects in Global for objFun -- this is only when not using multi-machine cluster
  advectionDir <- params$advectionDir
  advectionMag <- params$advectionMag
  omitPastPines <- FALSE # was TRUE ... but maybe bad because there are many that have multiple years
  dispersalKernel <- "Weibull"
  # dispersalKernel <- "Exponential"
  p <- do.call(c, params[c("meanDist", "advectionMag", "sdDist")])#, "advectionDir")])
  # p["meanDist"] <- 1e4 # Average distance per year -- in m
  # p["advectionMag"] <- 5e3 # Average distance of the asymmetry -- in m
  ##p["advectionDir"] <- 90
  sdDist <- p["sdDist"] # sqrt(variance(of the Weibull distribution)) --> contributes to shape and scale parameters
  #p["meanDistSD"] <- 20 # Interannual variation # not needed once we have wind and climate suitability
  #p["advectionDirSD"] <- 20
  # lower <- c(25000,  50,  0.9, 1.001)
  # upper <- c(65000, 2000, 2.5, 1.7)
  lower <- c(5000,  0.1,  0.5)
  upper <- c(105000, 2000, 2.5)
  p[] <- sapply(seq_along(p), function(x) runif(1, lower[x], upper[x]))

  fitType <- "distanceFromEachPoint"
  subsampleFrom <- NULL#3
  libPaths <- .libPaths()

  # Identify the objsToExport -- needs to be ALL the formalArgs of the objFun
  #   however, don't pass ones that are passed manually at the DEoptim call
  #   Generally, only pass "very small" objects in the DEoptim call. Large ones do this way.
  objsToExport <- setdiff(formalArgs("objFun"),
                          # The following are small and passed at the time of the DEoptim call
                          c("p", "reps", "quotedSpread", "fitType", "distanceFunction"))
  # need both the objsToExport vector, and the object called "objsToExport"
  objsToExport <- c(objsToExport, "reqdPkgs", "objsToExport", "libPaths")

  # put objsToExport into .GlobalEnv -- the DEoptim function gets all arguments from .GlobalEnv
  list2env(mget(objsToExport), envir = .GlobalEnv)

  # 213 appears slower than rest
  ipsNum <- c(189, 184, 220, 213, 217, 97, 106)
  ips <- paste0("10.20.0.", ipsNum)
  # ips <- c("localhost", ips)
  #ips <- c("10.20.0.184", "localhost", "10.20.0.220", "10.20.0.213", "10.20.0.217", "10.20.0.97") ## TODO: don't hardcode this. pass as param
  #ips <-   rep(ips, c(29, 16, 20, 0, 16, 17, 13))
  # ips <- rep(ips, c(16, 27, 25, 6, 16, 10))
  #                c(19, 22, 24, 6, 19, 10)
  #ips <- rep(ips, c(27, 16, 22, 19, 11, 5))
  #print(length(ips))
  #adjustments <- c(1.25, 1.1, #0.8,
  #                 0.8, 1.4) # relative, manual, highly idiosyncratic, depends on current use of machines
  #ips <- c("10.20.0.106", "10.20.0.68", "10.20.0.213")
  adjustments = 1#c(0.9, 2, 0.7, 2, 0.4, 0.4)
  numCoresNeeded <- 110 # This is one per true core on 4 machines, 56, 56, 28, 28 cores each

  # ips <- "10.20.0.213"
  # adjustments = 1
  # numCoresNeeded <- 3 # This is one per true core on 4 machines, 56, 56, 28, 28 cores each
  reqdPkgs <- grep(paste(collapse = "|", c("SpaDES.tools", "raster", "CircStats", "data.table", "purrr")), reqdPkgs, value = TRUE)
  if (isTRUE(type == "DEoptim")) {
    cl <- clusterSetup(workers = ips, objsToExport = objsToExport,
                       reqdPkgs = reqdPkgs, libPaths = libPaths, doSpeedTest = 1, # fn = fn,
                       # quotedExtra = quote(install.packages(c("rgdal", "rgeos", "sf",
                       #                                        "sp", "raster", "terra", "lwgeom"),
                       #                                      repos = "https://cran.rstudio.com")),
                       numCoresNeeded = numCoresNeeded, adjustments = adjustments)
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
                                               itermax = 50,
                                               NP = numCoresNeeded),
                     fitType = fitType,
                     subsampleFrom = subsampleFrom)

    saveRDS(DEout, file = file.path("outputs", paste0("DEout_", format(stPre), ".rds")))
    sim2 <- get("sim", whereInStack("sim"))
    parms <- params(sim2)
    saveRDS(parms, file = file.path("outputs", paste0("parms_", format(stPre), ".rds")))
    stPost <- Sys.time()
    print(difftime(stPost, stPre))
    try(parallel::stopCluster(cl))
    browser()
    sim$DEoptim <- DEoptim

  } else if (isTRUE(type == "runOnce")) {
    #profvis::profvis(
    st1 <- system.time({
      out <- objFun(quotedSpread = quotedSpread, reps = 1, p = p, fitType = fitType,
                    massAttacksStack = massAttacksStack, massAttacksDT = massAttacksDT, subsampleFrom = subsampleFrom)
    })
    #)
    print(st1)
  } else if (isTRUE(type == "optim")) {
    numCoresNeeded <- 5
    ips <- "10.20.0.213"
    cl <- clusterSetup(workers = ips, objsToExport = objsToExport,
                                reqdPkgs = reqdPkgs, libPaths = libPaths, doSpeedTest = 0, # fn = fn,
                                numCoresNeeded = numCoresNeeded)
    on.exit(try(parallel::stopCluster(cl)), add = TRUE)
    # cl <- parallel::makeForkCluster(5)
    stPre <- Sys.time()
    # optimOut <- optim(fn = objFun, par = p,
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
    try(parallel::stopCluster(cl))
    saveRDS(optimOut, file = file.path("outputs", paste0("optimOut_", format(stPre), ".rds")))
    stPost <- Sys.time()

  }
  objs <- dir("outputs", pattern = "DEout", full.names = TRUE)
  DEout <- readRDS(tail(objs, 1))
  pBest <- p
  pBest[] <- DEout$optim$bestmem[] # use [] to keep names


  return(invisible(pBest))
  #########################################################

  # visualizeKernel(p, massAttacksStack, maxDistance = maxDistance)
  ###############################################
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

# clusterSetupSingles <- function(workers, objsToExport, reqdPkgs,
#                          libPaths = .libPaths()[1], doSpeedTest = FALSE, envir = parent.frame(),
#                          fn = ".allObjs.rda", numCorlesNeeded, adjustments = rep(1, length(workers))) {
#
#   message("Starting cluster with 1 core per machine -- install reqdPkgs; copy objects; write to disk")
#   uniqueWorkers <- unique(workers)
#   if (identical(workers, uniqueWorkers) && length(workers) < numCoresNeeded) {
#
#     clSingle <- future::makeClusterPSOCK(workers = uniqueWorkers, revtunnel = TRUE)
#     on.exit(try(parallel::stopCluster(clSingle), silent = TRUE), add = TRUE)
#     # NumPopulations <- 118
#
#     out2 <- determineClusters(clSingle = clSingle, doSpeedTest = doSpeedTest,
#                               uniqueWorkers = uniqueWorkers, numCoresNeeded = numCoresNeeded,
#                               adjustments = adjustments)
#     reproducible::messageDF(out2)
#
#     clusterIPs <- rep(out2$.id, out2$cores)
#   } else {
#     clusterIPs <- if (length(workers) > numCoresNeeded)
#       sample(workers, size = numCoresNeeded)
#     else
#       workers
#   }
#   # parallel::stopCluster(clSingle)
#
#   # clSingle <- future::makeClusterPSOCK(workers = unique(clusterIPs), revtunnel = TRUE)
#   clusterExport(clSingle, varlist = objsToExport, envir = envir)
#   clusterExport(clSingle, varlist = c("fn", "reqdPkgs", "libPaths", "objsToExport"), envir = environment())
#
#   clusterEvalQ(clSingle, {
#     if (any(!dir.exists(libPaths)))
#       dir.create(libPaths[1], recursive = TRUE)
#     .libPaths(libPaths)
#     ## set CRAN repos; use binary linux packages if on Ubuntu
#     local({
#       options(Ncpus = parallel::detectCores() / 2)
#       options("repos" = c(CRAN = "https://cran.rstudio.com"))
#
#       if (Sys.info()["sysname"] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
#         .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
#         .user.agent <- paste0(
#           "R/", getRversion(), " R (",
#           paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"]),
#           ")"
#         )
#         options(repos = c(CRAN = paste0("https://packagemanager.rstudio.com/all/__linux__/",
#                                         .os.version, "/latest")))
#         options(HTTPUserAgent = .user.agent)
#       }
#     })
#
#     if (!suppressWarnings(require("Require"))) {
#       install.packages("Require")
#     }
#     # install.packages(c("rgdal", "rgeos", "sf", "sp", "raster", "terra", "lwgeom"), repos = "https://cran.rstudio.com")
#     suppressMessages(Require::Require(reqdPkgs, install = TRUE, require = FALSE))
#     save(list = objsToExport, file = fn)
#   })
#   parallel::stopCluster(clSingle)
#
#   return(clusterIPs)
# }

# determineClusters <- function(clSingle, doSpeedTest, uniqueWorkers, numCoresNeeded, adjustments) {
#
#   if (isTRUE(doSpeedTest) && length(unique(workers)) > 1) {
#     message("testing speed on each to estimate number cores to use")
#     out <- clusterEvalQ(clSingle, {
#       ss <- system.time({
#         for (i in 1:10000) rnorm(1e4)
#       })
#       data.table::data.table(elapsed = ss[3], trueCores = parallel::detectCores())
#     })
#   } else {
#     out <- clusterEvalQ(clSingle, {
#       data.table::data.table(elapsed = 1, trueCores = parallel::detectCores())
#     })
#   }
#
#   # parallel::clusterEvalQ(clSingle, system("pkill -f workRSOCK"))
#   names(out) <- uniqueWorkers
#   out2 <- rbindlist(out, idcol = TRUE)
#   relSpeed <- out2$trueCores/out2$elapsed*numCoresNeeded
#   relSpeed <- relSpeed/numCoresNeeded
#   relSpeed <- relSpeed/max(relSpeed)
#   out2[, relSpeed := relSpeed]
#   nonHTcores <- out2$trueCores/2
#   out2[, nonHTcores := nonHTcores]
#   sumNonHTcores <- sum(nonHTcores)
#   # needHTcores <- max(numCoresNeeded, numCoresNeeded - sumNonHTcores)
#
#   out2[, cores := round(nonHTcores * adjustments / max(adjustments))]
#
#   # Increase ncores upwards
#   m <- 0
#   while (sum(out2$cores) < numCoresNeeded) {
#     m <- ((m + 1) - 1) %% NROW(out2) + 1
#     out2[m, cores := cores + relSpeed]
#   }
#   # Decrease down to 1 each
#   out2[, cores := ceiling(cores)]
#   m <- 0
#   set(out2, NULL, "cores", as.numeric(out2$cores))
#   while ((sum(floor(out2$cores)) > numCoresNeeded) && any(out2$cores > 1)) {
#     m <- ((m + 1) - 1) %% NROW(out2) + 1
#     if (out2[m,]$cores > 1) {
#       maxC <- max(out2$cores)
#       out2[m, cores := cores - 1]
#     }
#   }
#   # set(out2, NULL, "cores", floor(out2$cores))
#   # Decrease down to 0 or 1 each
#   m <- 0
#   while (sum(out2$cores) > numCoresNeeded && any(out2$cores > 0)) {
#     m <- ((m + 1) - 1) %% NROW(out2) + 1
#     if (out2[m,]$cores > 0) {
#       out2[m, cores := cores - 1]
#     }
#   }
#
#   return(out2)
# }

grow <- function(sim) {
  ## determine the actual growth based on the actual number of attacked trees/ha
  sim$massAttacksDT[, greenTreesYr_t := xt(ATKTREES, CLIMATE, P(sim)$dataset,
                                           sim$massAttacksStack, mod$growthData)]

  return(invisible(sim))
}


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

# clusterSetup <- function(workers, objsToExport, reqdPkgs,
#                          libPaths = .libPaths()[1],
#                          doSpeedTest = FALSE, envir = parent.frame(),
#                          fn = ".allObjs.rda", numCoresNeeded = ceiling(detectCores() * 0.8),
#                          adjustments = rep(1, length(workers))) {
#   clusterIPs <- clusterSetupSingles(workers = workers, objsToExport = objsToExport,
#                              reqdPkgs = reqdPkgs, fn = fn, libPaths = libPaths, doSpeedTest = doSpeedTest,
#                              numCoresNeeded = numCoresNeeded)#, adjustments = adjustments)
#
#   message("Starting cluster with all cores per machine")
#   cl <- future::makeClusterPSOCK(clusterIPs, revtunnel = TRUE)
#   # cl <- future::makeClusterPSOCK(workers = clusterIPs, revtunnel = TRUE)
#   # on.exit(try(parallel::stopCluster(cl)), add = TRUE)
#   clusterExport(cl, varlist = c("fn", "reqdPkgs"), envir = environment())
#   st <- system.time(clusterEvalQ(cl, {
#     load(file = fn, envir = .GlobalEnv)
#     .libPaths(libPaths)
#     suppressMessages(
#       Require::Require(reqdPkgs, install = FALSE)
#     )
#   }))
#   st <- system.time(clusterEvalQ(cl, {
#     try(unlink(fn), silent = TRUE)
#   }))
#
#   return(cl)
# }

plotKernels <- function(p, reps = 10) {
  if ( (is.na(p["meanDistSD"])))
    reps = 1
  nrows <- max(1, floor(sqrt(reps + 1)))
  ncols <- max(1, ceiling((reps + 1)/nrows))
  par(mfrow = c(nrows, ncols))
  if (is.null(names(p))) stop("p must be named and have these 3: meanDist, sdDist, meanDistSD")
  for (i in 1:reps) {
    meanDist <- p["meanDist"]
    if (!is.na(p["meanDistSD"]))
      meanDist <- rlnorm(1, log(p["meanDist"]), log(p["meanDistSD"]))

    sdDist = p["sdDist"]
    shapeScale <- weibullShapeScale(meanDist = meanDist, sdDist = sdDist)
    # mn <- (meanDist)
    # sd <- mn/sdDist # 0.8 to 2.0 range
    # shape <- (sd/mn)^(-1.086)
    # scale <- mn/exp(lgamma(1 + 1/shape))

    ww <- rweibull(1e4, shape = shapeScale$shape, scale = shapeScale$scale)
    hist(ww, xlim = c(0, 14e4), main = "", xlab = "Distance (m) from source")
    abline(v = meanDist, col = "red")
    abline(v = p[["meanDist"]], col = "green")
  }
  plot(type = "n", x = 1:10, axes = F, xlab = "", ylab = "")
  legend(x = 2, y = 8, col = c("red", "green"), lty= 1,
         legend = c("annual average", "global average"),
         cex = 1.5)
  mtext(text = "Dispersal kernels", outer = TRUE, line = -3, xpd = TRUE)
}

kernelFn <- function(meanDist, sdDist, dispersalKernel, dist) {
  if (grepl("Weibull", dispersalKernel)) {
    shapeScale <- weibullShapeScale(meanDist, sdDist)

    if (grepl("Weibull", dispersalKernel)) {
      prob <- dweibull(dist, scale = shapeScale$scale, shape = shapeScale$shape)
    }

    infs <- is.infinite(prob)
    if (any(infs))
      prob[infs] <- 0
  } else {
    prob <- dexp(dist, rate = 1/meanDist)

  }
  prob
}

weibullShapeScale <- function(meanDist, sdDist) {
  mn <- (meanDist)
  sd <- mn/sdDist # 0.8 to 2.0 range
  shape <- (sd/mn)^(-1.086)
  scale <- mn/exp(lgamma(1+1/shape))
  return(list(scale = scale, shape = shape))
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
  raster::stack(stk)
}


stacksPredVObs <- function(massAttacksStack, predictedStack, propPineMap, threshold = 10, plot.it = FALSE) {
  yrNames <- names(massAttacksStack)
  yrNamesPlus1 <- yrNamesPlus1(yrNames)
  stk <- raster::stack()
  ROCs <- list()
  if (isTRUE(plot.it)) {
    numR <- floor(sqrt(length(yrNames)))
    par(mfrow = c(numR, ceiling(length(yrNames)/numR)))
  }
  for (lay in seq(yrNames)) {

    preds <- predictedStack[[yrNamesPlus1[lay]]][]
    # preds[preds < threshold] <- 0
    whHasVal <- which(!is.na(preds) |
            !is.na(massAttacksStack[[lay]][]) |
            propPineMap[] > 0)
    preds[is.na(preds)] <- 0

    # setColors(predictedStack[[yrNamesPlus1[lay]]], 8) <- "Reds"
    Both <- raster(massAttacksStack[[lay]])
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
        print(paste("year: ", yrNamesPlus1[lay],", threshold: ", threshold, ", ROC: ", round(ROCs[[lay]]$ROC, 2)))
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

plotCentroidShift <- function(centroids, title) {
  ggplot(centroids[], aes(x = Year, y = value, group = xy, color = type)) +
    geom_point() +
    facet_grid(vars(xy), scales = "free_y") +
    ggplot2::theme_bw() +
    ggtitle(title)
}

plotHistsOfDEoptimPars <- function(DEout, title) {
  (paramHists <- DEout %>%
    ggplot( aes(x=value, fill=Parameter)) +
    geom_histogram() +# color="#e9ecef", alpha=0.6, position = 'identity') +
    facet_grid(. ~ Parameter, scales = "free") +
    theme_bw() +
    ggtitle(title))
}
