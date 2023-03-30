crsLeaflet <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

if (FALSE) {
  saveSimList(sim, filename = "theSimList.qs")
  qs::qsave(sim$pineMap, "pineMap.qs")
  qs::qsave(sim$plotStack, "plotStack.qs")
  sim <- list()
  sim$plotStack <- qs::qread(file = "plotStack.qs")
  sim$pineMap <- qs::qread(file = "pineMap.qs")
  googledrive::drive_upload("pineMap.qs", path = "Shared/Hosted/leaflet/")
  googledrive::drive_upload("plotStack.qs", path = "Shared/Hosted/leaflet/")
  plotStackLeaf <- raster::stack(
    lapply(raster::unstack(sim$plotStack), function(r) {
      r <- projectInputs(r, targetCRS = crsLeaflet,
                  method = "ngb")
      r1 <- raster(r)
      r1[] <- r[]
      r1 <- writeOutputs(r1, filename2 = paste0(names(r), ".tif"))
    }))
  pineMap <- writeRaster(sim$pineMap, filename = "PineMap.tif", overwrite = TRUE)
}

# Require::Require(c("shiny", "leaflet", "leafem"))
#' RasterStacks are assumed to be annual time series, where the name of each layer
#' must have
#' @importFrom shiny server ui bootstrapPage absolutePanel colorNumeric
#' @importFrom raster minValue maxValue
#' @importFrom leaflet leafletOutput
#' @importFrom leafem addGeotiff
#' @param ... Named spatial objects e.g., \code{MPB = mpbStack}
#' @param opts a list of lists of options for each named object in \code{...}.
#'   These are passed to \code{addLegend} and \code{addGeotiff} and \code{addPolygons}
#'
#' @examples
opts = list(MPB = list(labels = c("LastYr", "NextYr", "LastAndNext", "PredYr", "LastYrAndPredYr",
                                  "NextYrAndPredYr", "LYRandNYandPY"),
                       colors = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "red", "red", "#E5C494")),
            PineMap = list(title = "Pine % cover"))
#' leafletSpaceTime(studyArea = sim$studyArea, MPB = plotStackLeaf, opts = opts)
leafletSpaceTime <- function(..., opts) {
  yrs <- 2011:2021
  ui <- bootstrapPage(
    tags$style(type = "text/css", "html, body {width:100%;height:100%}"),
    leafletOutput("mymap", width = "100%", height = "100%"),
    absolutePanel(top = 10, right = 10,
                  sliderInput("year", "Year", min(yrs), max(yrs),
                              value = yrs[1], step = 1, sep = "",
                  )#,
                  #selectInput("colors", "Color Scheme",
                  #            rownames(subset(brewer.pal.info, category %in% c("seq", "div")))
                  #),
                  #checkboxInput("legend", "Show legend", TRUE)
    )
  )

  # palGrey <- colorNumeric("Greys", 1:10*10)(10:100)

  objs <- list(...)

  classes <- list()
  isBrick <- which(sapply(objs, function(x) is(x, "RasterBrick")))
  if (any(isBrick)) stop("Cannot work with RasterBricks; please convert to RasterStack")
  classes$stack <- which(sapply(objs, function(x) is(x, "RasterStack")))
  classes$raster <- which(sapply(objs, function(x) is(x, "RasterLayer")))
  classes$vector <- which(sapply(objs, function(x)
    is(x, "sf") || is(x, "Spatial")))

  namesObjs <- names(objs)
  objs <- Cache(makeLongLat, ..., classes = classes)
  names(objs) <- namesObjs

  ranges <- lapply(objs, function(x) 0:1)
  colrs <- lapply(objs, function(x) "Reds")
  rasts <- c(classes$stack, classes$raster)
  if (length(rasts)) {
    ranges[rasts] <- lapply(rasts, function(clas)
      c(floor(min(minValue(objs[[clas]]))), ceiling(max(maxValue(objs[[clas]])))))
    colrs[rasts] <- lapply(names(objs[rasts]), function(x) colorNumeric("Reds", ranges[[x]] + c(-1, +1)))
  }

  colrsStatic <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "red", "red", "#E5C494")
  yrs <- format(Sys.Date(), format = "%Y")
  if (length(classes$stack))
    lapply(objs[classes$stack][1], function(s) {
      yrs <<- as.integer(gsub("[[:alpha:]]*", "", names(s)))
    })

  server <- function(input, output, session) {
    firstTime <- TRUE
    output$mymap <- renderLeaflet({
      ll <- leaflet() %>%
        addTiles()#
      if (any(classes$raster)) {
        out <- lapply(objs[classes$raster], function(rr) {
          # opts <- options(reproducible.cacheSpeed = "fast")
          browser()
          nBreaks <- 10
          brks <- c(0, seq(minValue(rr), maxValue(rr), nBreaks))
          palGrey <- rev(grey(scales::rescale(brks, c(2,10)/10)))
          palGrey[1] <- paste0(palGrey[1], "00")
          # brks <- seq(minValue(rr), maxValue(rr), nBreaks)
          ll <<- ll %>% addGeotiff(file = Filenames(rr), #x = rr, #
                                     opacity = 0.4, project = FALSE,
                       group = "Pine", # omitArgs = "map",
                       autozoom = if (isTRUE(firstTime)) {
                         firstTime <<- FALSE
                         TRUE
                       } else { FALSE },
                       colorOptions = leafem:::colorOptions(
                         palette = palGrey#hcl.colors(nBreaks, "inferno")#viridisLite::inferno,#colorRampPalette(grey(1,10/10)),#viridisLite::mako,#colorNumeric("Blue", domain = c(minValue(rr), maxValue(rr)))# viridisLite::inferno
                         , breaks = brks
                         , na.color = "transparent"
                       )
          ) %>%
            addLegend(position = "bottomright"
                      #, values = 1:10*10
                      , colors = palGrey#grey(1:10/10)
                      # , pal = viridisLite::inferno
                      , opacity = 1
                      , group = "Pine"
                      , labels = brks
                      , title = "Pine % cover")
          # options(opts)
        })
      }

      if (any(classes$vector)) {
        out <- lapply(objs[classes$vector], function(rr) {
          ll <<- ll %>% addPolygons(data = rr, fillColor = "transparent", weight = 1)
        })
      }

      ll <- ll %>%
        addLayersControl(
          # baseGroups = c("OSM (default)", "Toner", "Toner Lite"),
          overlayGroups = c("MPB", "Pine")
          , position = "topleft"
          # baseGroups = names(plotStack),
          , options = layersControlOptions(collapsed = FALSE)
        )

      #  addPolygons(fillColor = "transparent", weight = 1) %>%
      # ll %>% addLegend(position = "bottomright"
      #                  , colors = colrs
      #                  , opacity = 1
      #                  , labels = 1:7#raster::levels(plotStack[[1]])[[1]]$Label
      #                  , title = "Predictions vs Observations") # %>%
      # addLegend(position = "bottomright"
      #           , values = 1:10*10
      #           # , colors = grey(1:10/10)
      #           , pal = palGrey
      #           , opacity = 1
      #           # , labels = raster::levels(plotStack[[1]])[[1]]$Label
      #           , title = "Pine percent cover")

      if (any(classes$stack)) {
        out <- lapply(classes$stack, function(rrInd) {
          nam <- names(objs)[rrInd]
          optsHere <- modifyList(list(position = "bottomright"
                                      , colors = "Blues"
                                      , opacity = 1
                                      # , pal = colrs[[rrInd]]
                                      # , values = ranges[[rrInd]]
                                      # , labels = raster::levels(rr[[1]])[[1]]$Label
                                      , title = nam
                                      ),
                                 opts[[nam]])
          ll <<- SpaDES.tools:::docall(addLegend, append(list(map = ll), optsHere))
          # ll <<- ll %>% addLegend(position = "bottomright"
          #                         , colors = colrsStatic#"Blues"
          #                         , opacity = 1
          #                         # , pal = colrs[[rrInd]]
          #                         # , values = ranges[[rrInd]]
          #                         , labels = raster::levels(rr[[1]])[[1]]$Label
          #                         , title = "Predictions vs Observations")
        })
      }

      nam <- names(objs)[classes$stack]

      colOpts <- colorOptions(
        # palette = colrsStatic,
        na.color = "transparent"
      )
      colsOptsHere <- modifyList(colOpts, list(palette = opts[[nam]]$colors))

      rr <- objs$MPB$predX2011_lf
      # optsHere <- modifyList(list(#x = rr,
      #   file = Filenames(rr),#, omitArgs = c("map", "data"))
      #   autozoom = if (isTRUE(firstTime)) {
      #     firstTime <<- FALSE
      #     TRUE
      #   } else { FALSE }
      #   # , x = rr
      #   , layerId = paste0("ras", input$year)
      #   , opacity = 1
      #   , project = FALSE
      #   , colorOptions = colsOptsHere
      #   #, colorOptions = colorOptions(
      #   #  palette = colrsStatic,
      #   #  na.color = "transparent"
      #   #)
      #   , group = "MPB",
      #   title = nam
      # ),
      # opts[[nam]])
      # # lp <- leafletProxy("mymap", data = rr) #%>%
      # ll <<- SpaDES.tools:::docall(addGeotiff, append(list(map = ll), optsHere)) #%>%
        #clearGroup(nam)
      browser()
      ll <- ll %>% addGeotiff(#x = rr,#
        file = Filenames(rr),#, omitArgs = c("map", "data"))
        # leafem::addGeoRaster(map = lp, # omitArgs = c("map"),
        # leafem::addGeoRaster(#omitArgs = c("map"),
        autozoom = if (isTRUE(firstTime)) {
          firstTime <<- FALSE
          TRUE
        } else { FALSE }
        # , x = rr
        , layerId = paste0("ras", input$year)
        , opacity = 1
        , project = FALSE
        , colorOptions = colorOptions(
          domain = c(1,7),
          palette = colrsStatic,
          na.color = "white"
        )
        , group = "MPB")


      ll
    }) # end renderLeaflet

    YrRas <- reactive({
      if (length(classes$stack)) {
        yr <- as.integer(gsub("[[:alpha:]]*", "", names(objs[[classes$stack]])))
        lay <- grep(input$year, names(objs[[classes$stack]]), value = TRUE)
        out <- objs[[classes$stack]][[lay]] # TODO need layer names to match stack
        return(out)
      }
    })


    ### Stacks #################################################

    # observe({
    #
    #   rr <- YrRas()
    #   if (!is.null(rr) && FALSE) {
    #     # rr[555] <- 6
    #     browser()
    #     nam <- names(objs)[classes$stack]
    #
    #     colOpts <- colorOptions(
    #       # palette = colrsStatic,
    #       na.color = "transparent"
    #     )
    #     colsOptsHere <- modifyList(colOpts, list(palette = opts[[nam]]$colors))
    #
    #     optsHere <- modifyList(list(#x = rr,
    #                                 file = Filenames(rr),#, omitArgs = c("map", "data"))
    #                                 autozoom = if (isTRUE(firstTime)) {
    #                                   firstTime <<- FALSE
    #                                   TRUE
    #                                 } else { FALSE }
    #                                 # , x = rr
    #                                 , layerId = paste0("ras", input$year)
    #                                 , opacity = 1
    #                                 , project = FALSE
    #                                 , colorOptions = colsOptsHere
    #                                 #, colorOptions = colorOptions(
    #                                 #  palette = colrsStatic,
    #                                 #  na.color = "transparent"
    #                                 #)
    #                                 , group = "MPB",
    #                                 title = nam
    #     ),
    #     opts[[nam]])
    #     lp <- leafletProxy("mymap", data = rr) #%>%
    #     lp <<- SpaDES.tools:::docall(addGeotiff, append(list(map = lp), optsHere)) %>%
    #       clearGroup(nam)
    #     # lp <<- SpaDES.tools:::docall(addGeoRaster, append(list(map = lp), optsHere)) %>%
    #     #   clearGroup(nam)
    #
    #     # addGeotiff(map = lp, file = Filenames(rr),#, omitArgs = c("map", "data"))
    #     #            # leafem::addGeoRaster(map = lp, # omitArgs = c("map"),
    #     #            # leafem::addGeoRaster(#omitArgs = c("map"),
    #     #            autozoom = if (isTRUE(firstTime)) {
    #     #              firstTime <<- FALSE
    #     #              TRUE
    #     #            } else { FALSE }
    #     #            # , x = rr
    #     #            , layerId = paste0("ras", input$year)
    #     #            , opacity = 1
    #     #            , project = FALSE
    #     #            , colorOptions = colorOptions(
    #     #              palette = colrsStatic,
    #     #              na.color = "transparent"
    #     #            )
    #     #            , group = "MPB"
    #     # ) %>% clearGroup("MPB")
    #   }
    #   #browser()
    #   #lp %>% lp1 %>% clearGroup("MPB")
    #
    #
    #   # lp <- leafletProxy("mymap", data = YrRas())
    #   #
    #   # lapply(classes$stack, function(x) {
    #   #   opts <- options(reproducible.cacheSpeed = "fast")
    #   #   lp <<- Cache(leafem::addGeoRaster, map = lp, omitArgs = c("map", "autozoom")
    #   #                , autozoom = FALSE
    #   #                # if (isTRUE(firstTime)) {
    #   #                #   firstTime <<- FALSE
    #   #                #   TRUE
    #   #                # } else { FALSE }
    #   #                , x = YrRas()
    #   #                , layerId = paste0("ras", input$year)
    #   #                , opacity = 1
    #   #                , project = FALSE
    #   #                , colorOptions = colorOptions(
    #   #                  palette = colrsStatic#colrs[[x]]
    #   #                  , na.color = "transparent"
    #   #                )
    #   #                , group = "MPB"
    #   #   ) %>%
    #   #     clearGroup("MPB")
    #   #   options(opts)
    #   # })
    #   # # browser()
    #   # lp
    # })
  }

  shinyApp(ui, server)
}

makeLongLat <- function(..., classes) {
  objs <- list(...)

  objs[classes$vector] <- lapply(objs[classes$vector], function(p) {
    if (!is(p, "sf"))
      p <- sf::st_as_sf(p)
    if (!identical(st_crs(p), st_crs(4326)))
      p <- sf::st_transform( p, crs = 4326)
    return(p)
  })

  objs[classes$raster] <- lapply(objs[classes$raster], function(p) {
    fn <- reproducible:::.suffix(Filenames(p), suffix = "_lf")
    if (!identical(st_crs(p), st_crs(4326))) {
      p <- Cache(projectInputs, p, targetCRS = crsLeaflet)
      p <- writeOutputs(p, filename = fn, overwrite = TRUE)
    }
    return(p)
  })

  rasLFTmplate <- raster(extent(0, 1, 0, 1), res = 1, val = 1, crs = crsLeaflet)
  objs[classes$stack] <- lapply(objs[classes$stack], function(stk) {
    stkList <- lapply(raster::unstack(stk), function(p) {
      fn <- reproducible:::.suffix(Filenames(p), suffix = "_lf")
      if (inherits(p, "Raster")) {
        p = stars::st_as_stars(p)
      }
      if (!sf::st_is_longlat(p)) {
        p = stars::st_warp(p, crs = 4326)
      }
      if (inherits(p, "stars_proxy")) {
        file.copy(p[[1]], fn)
      }
      if (!inherits(p, "stars_proxy")) {
        stars::write_stars(p, dsn = fn)
      }
      p <- raster::raster(fn)
      p <- setMinMax(p)
      return(p)
    })
    raster::stack(stkList)

  })

  objs
}

if (FALSE) {
  leafletSpaceTime(prediction = plotStackLeaf, studyArea = sim$studyArea, pineMap = sim$pineMap)
}
