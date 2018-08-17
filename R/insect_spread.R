#install_github("PredictiveEcology/SpaDES.tools@prepInputs")
#devtools::load_all("~/GitHub/SpaDES.tools")



#' Simulate asymmetric propagule spread
#'
#' Propagules distributed based on 1) distance from focal cell (following a dispersal
#' kernel -- currently negative exponentiol); 2) direction; 3) resistance surface.
#' Stop-based rule based on depletion of propagules from focal cells (loci).
#'
#' Produces ovoid dispersal "clouds", instead of cones. # TODO: come back to this
#'
#' @param r      A resistance surface (i.e., habitat quality) \code{RasterLayer} map,
#'               whose values range between 0 and 1, with larger values indicating
#'               increased likelihood of propagule settling (i.e., resistance to movement).
#'
#' @param loci   Integer vector of indices into \code{r} indicating start locations.
#'
#' @param asymmetry Ratio indicating the "strength" of the asymmetry to use for spread.
#'                  Larger numbers indicate greater asymmetry.
#'
#' @param asymmetryAngle Indicates the direction of spread bias. North is 0 degrees.
#'
#' @param lambda parameter in the negative exponential dispersal kernal function
#'               describing the steepness of curve.
#'               Greater values (towards 1) is steeper, lower (towards 0) is flatter.
#'
#' @param saturationDensity DESCRIPTION NEEDED
#'
#' @param total  The maximum number of propagules per pixel.
#'
#' @param Ncpus  Maximum number of processor cores to use.
#'
#' @param maxIterations The maximum number of iterations, regardless of all other drivers
#'
#' @param ... Any other arguments passed to spread2
#'
#' @param debug  Logical indicating whether debugging/status information should be
#'               printed to the console. Default \code{FALSE}.
#'
#' @return A \code{RasterLayer} whose values indicate the number of propagules in each pixel.
#'
#' @author Eliot McIntire
#' @export
#' @importFrom data.table := rbindlist set setDTthreads setkey
#' @importFrom raster raster
#' @importFrom SpaDES.tools spread2
#'
#' @example
#' library(quickPlot)
#' library(raster)
#' library(SpaDES.tools)
#'
#' # create a dummy raster
#' XMAX <- YMAX <- 400
#' r <- raster(extent(0, XMAX, 0, YMAX), res = 1)
#' pine <- gaussMap(r)
#' pine[] <- pine[] / maxValue(pine)
#'
#' r2 <- insect_spread(pine, loci = sample(ncell(r), 100),
#'                     asymmetry = 5, asymmetryAngle = 90,
#'                     lambda = 0.12,
#'                     saturationDensity = 332,
#'                     total = round(1125 * (250 / 100)^2),
#'                     debug = TRUE)
#'
#' # plot the result
#' dev.useRSGD(FALSE)
#' dev()
#' clearPlot()
#' Plot(r2, zero.color = "blue", new = TRUE, col = "Reds")
#'
insect_spread <- function(r, loci, asymmetry, asymmetryAngle, lambda,
                          saturationDensity, total,
                          Ncpus = getOption("Ncpus", parallel::detectCores() / 2),
                          debug = FALSE, maxIterations = 1e6, ...) {

  ## TODO: implement other dispersal kernels
  dispKern <- function(disFar, disNear, lambda) {
     (1 - exp(-lambda * disFar)) - (1 - exp(-lambda * disNear))
  }

  setDTthreads(Ncpus)

  pix <- r[]

  out <- spread2(r, start = loci, spreadProb = 1, asRaster = FALSE, iterations = 0,
                 asymmetry = asymmetry, asymmetryAngle = 0,
                 circle = TRUE,
                 returnDistances = TRUE, returnFrom = TRUE)
  set(out, , "abundanceActive", 0)
  set(out, , "abundanceSettled", 0)
  set(out, , "abundanceReceived", 0)
  set(out, , "Total", total)

  ## within a season, disperse the beetles
  done <- FALSE
  while (!done) {
    nrowOut <- NROW(out)
    out <- spread2(r, start = out, spreadProb = 1, asRaster = FALSE, iterations = 1,
                   asymmetry = asymmetry, asymmetryAngle = asymmetryAngle,
                   circle = TRUE, allowOverlap = NA,
                   returnDistances = TRUE, returnFrom = TRUE, ...)
    if (nrowOut != NROW(out)) {
      # might spread and have none that are <1 unit. Just go back to while loop
      set(out, , "order", seq_len(NROW(out)))
      attribs <- attr(out, "spreadState")

      setkey(out, initialPixels, from)
      out[, beetleSource := !is.na(abundanceActive),  by = c("initialPixels", "from")]
      receivers <- out[beetleSource == FALSE]
      sources <- out[beetleSource == TRUE]

      setkey(sources, initialPixels, pixels)
      setkey(receivers, initialPixels, from)
      outW_i_columns <- receivers[sources]
      outWLag1B <- outW_i_columns[beetleSource == FALSE]

      # Calculate the abundance received, as a function of distance
      outWLag1B[, abundanceReceived :=
                  ceiling(dispKern(effectiveDistance, i.effectiveDistance, lambda) *
                            proportion *
                            i.Total + i.abundanceActive *
                            proportion / sum(proportion)),
                by = c("initialPixels", "from")] # Extra ones if they didn't settle
      # Calculate the abundance received, as a function of angle,
      # which was already calculated in spread2, and is called "proportion".
      # The pmin is about saturation density.
      outWLag1B[, abundanceSettled := pmin(floor(abundanceReceived * pix[pixels] *
                                                   saturationDensity /
                                                   sum(abundanceReceived * pix[pixels])),
                                           ceiling(abundanceReceived * pix[pixels])),
                by = c("pixels")]
      outWLag1B[, abundanceActive := floor(abundanceReceived - abundanceSettled)]

      set(outWLag1B, , grep(colnames(outWLag1B), pattern = "^i\\.", value = TRUE), NULL)

      out <- rbindlist(list(sources, outWLag1B), fill = TRUE)

      # Squash down duplicates
      outSum <- out[, list(abundanceActive = sum(abundanceActive),
                           abundanceSettled = sum(abundanceSettled),
                           abundanceReceived = sum(abundanceReceived)),
                    by = c("initialPixels", "pixels")]
      out <- unique(out, by = c("initialPixels", "pixels"))
      set(out, , "abundanceActive", outSum$abundanceActive)
      set(out, , "abundanceSettled", outSum$abundanceSettled)
      set(out, , "abundanceReceived", outSum$abundanceReceived)
      set(out, , "Total", total)
      setattr(out, "spreadState", attribs)

      # Because of multiple starting loci, a given pixel can get overfull.
      overfull <- out[, sum(abundanceSettled), by = "pixels"]
      if (isTRUE(any(overfull$V1 > saturationDensity))) {
        pixelsWithTooMany <- overfull[V1 > saturationDensity]$pixels
        set(out, , "cumsumAbSett", 0)
        set(out, , "abundanceSettledTemp", 0)
        out[(pixels %in% pixelsWithTooMany), cumsumAbSett := cumsum(abundanceSettled),
            by = c("pixels")]
        cumsumAboveSaturation <- which(out$cumsumAbSett > saturationDensity)
        abundanceSettledTemp <- pmax(0, saturationDensity -
                                       (out$cumsumAbSett[cumsumAboveSaturation] -
                                          out$abundanceSettled[cumsumAboveSaturation]))
        set(out, cumsumAboveSaturation, "abundanceActive", out$abundanceActive[cumsumAboveSaturation] +
              (out$abundanceSettled[cumsumAboveSaturation] - abundanceSettledTemp))
        set(out, cumsumAboveSaturation, "abundanceSettled", abundanceSettledTemp)
      }

      if (all(out[, sum(abundanceSettled, na.rm = TRUE) >= unique(Total), by = "initialPixels"])) {
        done <- TRUE
      }

      if (attr(out, "spreadState")$totalIterations > maxIterations) {
        done <- TRUE
      }

      if (isTRUE(debug)) {
        ## print debugging/status info
        print(paste(attr(out, "spreadState")$totalIterations,
                    paste(out[, sum(abundanceSettled, na.rm = TRUE) >= unique(Total),
                              by = "initialPixels"]$V1,
                    collapse = ", ")))
      }
    }
  }
  aa <- raster(r)
  aa[] <- NA_integer_
  out1 <- out[, list(abundanceSettled = sum(abundanceSettled),
                     abundanceReceived = sum(abundanceReceived)), by = c("pixels")]
  aa[out1$pixels] <- out1$abundanceSettled

  return(aa)
}
