quotedSpreadOld <- quote({
  SpaDES.tools::spread3(start = starts,
                        rasQuality = propPineRas,
                        rasAbundance = massAttacksRas,
                        advectionDir = rnorm(1, p[3], p[6]),
                        advectionMag = p[2],
                        meanDist = rlnorm(1, log(p[[1]]), log(p[[5]])),
                        sdDist = p[4],
                        dispersalKernel = dispersalKernel,
                        plot.it = FALSE,
                        minNumAgents = minNumAgents,
                        verbose = 0,
                        skipChecks = TRUE,
                        saveStack = NULL)
})

# if (exists("EliotTesting")) {
#   EliotTesting <- TRUE
#   sim@params$mpbRedTopSpread$bgSettlingProp <- 0.7
#   sim@params$mpbRedTopSpread$advectionMag <- 20000
#   # minNumAgents <- 200
#   sim@params$mpbRedTopSpread$.plotInitialTime <- NA
# } else {
#   EliotTesting <- FALSE
# }
# propPineRas[] <- pmin(1, propPineRas[] + bgSettlingProp) ## TODO: why divide by 100??

# if (EliotTesting) { # TODO -- delete EliotTesting when no longer desired
#   a <- extent(studyArea)
#   starts <- massAttacksDT[ATKTREES > 0]$ID
#   d <- raster(massAttacksRas)
#   d[starts] <- 1
#   d <- crop(d, a)
#   starts <- which(d[] > 0)
#   massAttacksRas <- crop(massAttacksRas, a) * 300
#   propPineRas <- crop(propPineRas, a)
#   saveStack <- raster::rasterTmpFile()
# } else {
#   minNumAgents <- 50
#   starts <- massAttacksDT[["ID"]][massAttacksDT$ATKTREES > 0]
#   saveStack <- NULL
#   massAttacksRas <- massAttacksRas
# }

# mam <- raster::unstack(massAttacksStack)
# names(mam) <- names(massAttacksStack)
# massAttacksDTAllYears <- lapply(mam, function(x) {
#   pixels = which(x[] > 0)
#   setDT(list(pixels = pixels, ATKTREES = x[][pixels]))
# })
