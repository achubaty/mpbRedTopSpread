visualizeKernel <- function(p, massAttacksStack, maxDistance) {
  browser()





  maxD <- maxDistance + 125
  r <- raster(extent(-maxD+ 250, maxD + 250 , -maxD, maxD), res = 250)
  r[] <- 1
  r <- gaussMap(r)
  r <- r/maxValue(r)
  abund <- raster(r)
  #abund[mid] <- 10000

  spiralsD <- spiralDistances(r, maxDis = maxDistance, cellSize = res(massAttacksStack)[1])
  tos <- ( spiralsD[, c("row", "col")]) * res(r)
  colnames(tos) <- c("x", "y")

  p_sdDist = p[4]
  par(mfrow = c(4,4))
  for (i in 1:4) {
    mid <- SpaDES.tools::middlePixel(abund)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    p_meanDist = rlnorm(1, log(p[[1]]), log(p[[5]]))
    windDir <- rnorm(1, p["p_advectionDir"], p["advectionDirSD"]*3)
    print("aaa")
    dd <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                           p_meanDist = p_meanDist, p_sdDist = p["p_sdDist"], landscape = r,
                           propPineRas = r, fromCell = mid,
                           toCells = toCells,
                           p_advectionMag = p["p_advectionMag"], angle = ddd[, "angles"],
                           windDir = windDirRas[])
    abund[toCells] <- a1 <- -dd * 10000/(sum(-dd, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd1 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                            p_meanDist = p_meanDist, p_sdDist = p["p_sdDist"], landscape = r,
                            propPineRas = r, fromCell = mid,
                            toCells = toCells,
                            p_advectionMag = p["p_advectionMag"], angle = ddd[, "angles"],
                            windDir = windDirRas[])
    abund[toCells] <- a2 <- -dd1 * 10000/(sum(-dd1, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd2 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                            p_meanDist = p_meanDist, p_sdDist = p["p_sdDist"], landscape = r,
                            propPineRas = r, fromCell = mid,
                            toCells = toCells,
                            p_advectionMag = p["p_advectionMag"], angle = ddd[, "angles"],
                            windDir =  windDirRas[])
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
  maxD <- maxDistance + 125
  r <- raster(extent(-maxD+ 250, maxD + 250 , -maxD, maxD), res = 250)
  r[] <- 1
  r <- gaussMap(r)
  r <- r/maxValue(r)
  abund <- raster(r)
  #abund[mid] <- 10000

  spiralsD <- spiralDistances(r, maxDis = maxDistance, cellSize = res(massAttacksStack)[1])
  tos <- ( spiralsD[, c("row", "col")]) * res(r)
  colnames(tos) <- c("x", "y")

  p_sdDist = p[4]
  # p["p_advectionMag"] <- 15000
  par(mfrow = c(4,4))
  for (i in 1:4) {
    mid <- SpaDES.tools::middlePixel(abund)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    p_meanDist = rlnorm(1, log(p[[1]]), log(p[[5]]))
    windDir <- rnorm(1, p["p_advectionDir"], p["advectionDirSD"]*3)
    dd <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                           p_meanDist = p_meanDist, p_sdDist = p["p_sdDist"], landscape = r,
                           propPineRas = r, fromCell = mid,
                           toCells = toCells,
                           p_advectionMag = p["p_advectionMag"], angle = ddd[, "angles"],
                           windDir = windDirRas[])
    abund[toCells] <- a1 <- -dd * 10000/(sum(-dd, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd1 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                            p_meanDist = p_meanDist, p_sdDist = p["p_sdDist"], landscape = r,
                            propPineRas = r, fromCell = mid,
                            toCells = toCells,
                            p_advectionMag = p["p_advectionMag"], angle = ddd[, "angles"],
                            windDir = windDirRas[])
    abund[toCells] <- a2 <- -dd1 * 10000/(sum(-dd1, na.rm = TRUE))
    plot(abund)


    mid <- SpaDES.tools::middlePixel(abund) + sample(size = 1, -50000:50000)
    (midXY <- xyFromCell(r, mid))
    ddd <- distanceFromEachPoint(from = midXY, to = tos, angles = TRUE)
    nas <- is.na(toCells)
    ddd <- ddd[!nas,]
    toCells <- cellFromXY(r, xy = ddd[, c("x", "y")])
    dd2 <- distanceFunction(ddd[, "dists"], dispersalKernel = "exponential",
                            p_meanDist = p_meanDist, p_sdDist = p["p_sdDist"], landscape = r,
                            propPineRas = r, fromCell = mid,
                            toCells = toCells,
                            p_advectionMag = p["p_advectionMag"], angle = ddd[, "angles"],
                            windDir =  windDirRas[])
    abund[toCells] <-  a3 <- -dd2 * 10000/(sum(-dd2, na.rm = TRUE))
    plot(abund)
    a4 <- rescale(a1, c(0,1)) + rescale(a2, c(0,1)) + rescale(a3, c(0,1))
    abund[toCells] <- a4
    plot(abund)
  }
# Visualize with maps
# Kernel
advectionMagTmp <- p[2] ;
p_advectionDir <- p[3]
r <- raster(extent(-20000, 20000, -20000, 20000), res = 100)
r[] <- 1
abund <- raster(r)
mid <- SpaDES.tools::middlePixel(abund)
abund[mid] <- 10000
out <- spread3(mid, rasQuality = r, rasAbundance = abund,
               p_advectionDir = p_advectionDir,
               p_advectionMag = advectionMagTmp,
               p_meanDist = p_meanDist,
               p_sdDist = p_sdDist, dispersalKernel = "weibull", plot.it = F)
abund[out$pixels] <- out$abundSettled
abund[mid] <- max(abund[], na.rm = TRUE) * 1.5
clearPlot(); Plot(abund)


nams <- names(massAttacksStack)
# diffs <- list()
diffs <- lapply(seq(nams)[-length(nams)], function(i) {
  a <- massAttacksStack[[nams[i]]]
  b <- propPineRas
  whTos <- which(b[] >= 0.2)
  whFroms <- which(a[] > 0)
  if (length(whFroms) && length(whTos)) {
    env <- new.env()
    env$atksRasNextYr <- b
    env$atksKnownNextYr <- data.table(pixels = whTos)
    env$starts <- whFroms
    env$massAttacksRas <- a
    env$dispersalKernel <- "exponential"
    env$propPineRas <- propPineRas
    env$p <- p
    st <- system.time(out <- eval(quotedSpread, envir = env))
    pixels <- cellFromXY(massAttacksStack, out[, c("x", "y")])
    predictedMap <- raster(massAttacksStack)
    predictedMap[whTos] <- log(-out[, "val"] + 1)
    clearPlot()
    massAttacksRas <- env$massAttacksRas
    nextYearAttacks <- massAttacksStack[[nams[i+1]]]
    Plot(massAttacksRas, nextYearAttacks, predictedMap)


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

mn <- (p_meanDist + advectionMagTmp)
sd <- mn/p_sdDist # 0.8 to 2.0 range
shape <- (sd/mn)^(-1.086)
scale <- mn/exp(lgamma(1 + 1/shape))
aa <- curve(dweibull(x, shape = shape, scale = scale), from = 0, 1e5)

mn <- max(10, (p_meanDist - advectionMagTmp))
sd <- mn/p_sdDist # 0.8 to 2.0 range
shape <- (sd/mn)^(-1.086)
scale <- mn/exp(lgamma(1 + 1/shape))
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
  rr <- raster::raster(extent(-maxDist, maxDist, -maxDist, maxDist), res = res(propPineRas))
  out3 <- out2[abundSettled > 0]
  xx = out3$distance * (sin(out3$direction))
  yy = out3$distance * (-cos(out3$direction))
  pixels <- cellFromXY(rr, cbind(xx, yy))
  rr[pixels] <- out3$abundSettled
  raster::plot(rr)


  out2 <- out[, list(abundSettled = sum(abundSettled)), by = c("pixels")]
  rr <- raster::raster(propPineRas)
  rr[out2$pixels] <- out2$abundSettled
  rr2 <- trim(rr)
  rr3 <- raster(rr2)
  rr3[rr2[] > 10] <- 2
  rr4 <- trim(rr3)
  rr5 <- crop(rr2, rr4)
  Plot(rr5)


}
}
