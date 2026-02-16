
#' Show information in setA class \pkg{bamm}.
#' @importFrom methods new
#' @importFrom utils head
#' @param object An object of class setA
#' @importFrom methods show
#' @rdname show
#' @return Display information about the setA object
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export

methods::setMethod(f = "show",
                   signature = "setA",
                   function(object) {
                     slotsin <- methods::slotNames(object)

                     cat("Set A of the BAM digram it contains",
                         length(slotsin),"slots \n\n")
                     #cat("@niche_model, @suit_threshold,
                     #@cellIDs,@suit_values\n",
                     #   "@sparse_model,@coordinates,@eigen_vec,
                     #@eigen_val",sep = "")
                     npixs <- length(object@cellIDs)
                     if("niche_model" %in% slotsin){
                       cat("@niche_model: a niche model:\n\n")
                       print(object@niche_model)
                     }
                     if("suit_threshold" %in% slotsin){
                       cat("@suit_threshold: Threshold value used",
                           "to binarize model")
                       #print(object@suit_threshold)
                     }
                     if("cellIDs" %in% slotsin){
                       cat("@cellIDs: ids of the cells that have values",
                           paste0("(",npixs," pixels)"),"\n\n")
                     }
                     if("suit_values" %in% slotsin){
                       cat("@suit_values:",
                           "Suitability values of the continuous model\n\n")

                     }
                     if("sparse_model" %in% slotsin){
                       cat("@sparse_model:",
                           "A sparse square matrix of ",
                           npixs,"x",npixs,
                           "entries \n\n")

                     }
                     if("suit_values" %in% slotsin){
                       cat("@coordinates:",
                           "Pixel centroid coordinates of the model\n\n")

                     }



                    # if("occs_sparse" %in% slotsin){
                    #   cat("Number of occurences:",
                    #       (object@n_occs))
                    # }

                   })

#' Show information in csd class \pkg{bamm}.
#' @importFrom methods new
#' @param object An object of class setA
#' @rdname show
#' @return Display information about the csd object
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export

methods::setMethod(f = "show",
                   signature = "csd",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("Object of class csd it contains",
                         length(slotsin),"slots: \n\n")
                     cat("@connections: Geographic clusters data.frame \n\n")
                     print(head(object@connections,4))
                     cat("@interactive_map: A leaflet map showing the",
                         "geographic clusters\n\n")
                     cat("@raster_map: A raster map of the clusters")
                   })


#' Show information in pam class \pkg{bamm}.
#' @importFrom methods new
#' @param object An object of class pam
#' @return Display information about the pam object
#' @rdname show
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export

methods::setMethod(f = "show",
                   signature = "pam",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("Object of class pam it contains",
                         length(slotsin),"slots: \n\n")
                     cat("@pams: A list of sparse PAM matrices\n\n")
                     print(head(object@pams[[1]]))
                     cat("\n@sp_names: the species in the pam with ",
                         length(object@sp_names), "species:\n",
                         paste0(object@sp_names,collapse = ", "))
                     cat("\n\n@which_steps: time steps",object@which_steps,
                         "\n\n")
                     cat("@grid:", "raster grid of the studied area\n\n" )
                     print(object@grid)
                     cat("\n@cellIDs:", "site ids regarding the raster grid",
                         "\n\n" )

                   })


#' Show information in pam class \pkg{bamm}.
#' @importFrom methods new
#' @param object An object of class bioindiex_sparse
#' @rdname show
#' @return Display information about the bioindex_spars object
#' @export

methods::setMethod(f = "show",
                   signature = "bioindex_sparse",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("Object of class bioindex it contains",
                         length(slotsin),"slots: \n\n")
                     cat("@alpha: A sparse matrix with the richness of species",
                         "per site \n\n")
                     print(object@alpha,6)
                     cat("\n")
                     cat("@omega: A sparse matrix with the range size of every",
                         "species \n\n")
                     print(object@omega,6)
                     cat("\n")
                     cat("@dispersion_field: A sparse with the set of ranges",
                         "of all species that occur in at each locality \n\n")
                     print(object@omega,6)

                   })


#' Show information in setA class \pkg{bamm}.
#' @importFrom methods new
#' @param object An object of class setM
#' @return Display information about the setM object
#' @rdname show
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export

methods::setMethod(f = "show",
                   signature = "setM",
                   function(object) {
                     slotsin <- methods::slotNames(object)

                     cat("Set M of the BAM digram it contains",
                         length(slotsin),"slots \n\n")


                     cat("@coordinates: A matrix with longitude and latitude",
                         "values of each cell of the raster area\n\n")
                     print(head(object@coordinates))

                #cat("@initial_points: A list of inital coordinates where the",
                    #      "invasion process starts\n\n")
                     #if(length(object@initial_points)>0L)
                    #   print(head(object@initial_points))

                     cat("@eigen_val: Eigen values of the connectivity matrix",
                         "M\n\n")
                     if(length(object@eigen_val)>0L)
                       print(head(object@eigen_val))

                     cat("@eigen_vec: Eigen vector of the connectivity matrix",
                         "M\n\n")
                     if(length(object@eigen_vec)>0L)
                       print(head(object@eigen_vec))
                   })


#' Show information in diversity_range class \pkg{bamm}.
#' @importFrom methods new
#' @param object An object of class diversity_range
#' @rdname show
#' @return Display information about the diversity_range object
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export

methods::setMethod(f = "show",
                   signature = "diversity_range",
                   function(object) {
                     slotsin <- methods::slotNames(object)
                     cat("diversity_range object it contains",
                         length(slotsin),"slots \n\n")


                     cat("@alpha: A column vector of size ",object@nsites,
                         " with values of alpha diversity at each site",
                         "\n\n")
                     print(head(object@alpha))

                     cat("@alpha_raster: Alpha diversity raster",
                         "\n\n")
                     print(object@alpha_raster)

                     cat("@omega: A column vector of size ",object@nsps,
                         " with the range size of each species",
                         "\n\n")
                     print(head(object@omega))
                     cat("@dispersion_field: A column vector of size ",
                         object@nsites,
                         "with values of dispersion field at each site",
                         "\n\n")
                     print(head(object@dispersion_field))
                     cat("@dispersion_field_raster: Dispersion field raster",
                         "\n\n")
                     print(object@dispersion_field_raster)
                     cat("@null_dispersion_field_dist: null dispersion field",
                         "distribution. ")
                     cat("It is matrix of\n n =",object@nsites,"sites x","m =",
                         object@n_iterations, "simulations",
                         "used to generate random values of dispersion field" ,
                         "\n\n")
                     if(ncol(object@null_dispersion_field_dist)>=2){
                       print(head(object@null_dispersion_field_dist[c(1,2),
                                                                    c(1,2)]))
                     }

                     cat("@diversity_range_raster: Raster with diversity range",
                         "categories",
                         "\n\n")
                     print(object@diversity_range_raster)
                     cat("@xy_coordinates: Geographical coordinates of the",
                         "sites","\n\n")
                     print(head(object@xy_coordinates))

                   })


#if (!isGeneric("plot")) {setGeneric("plot", function(x,y,...)
#standardGeneric("plot"))}


#' Plot method for objects of class diversity_range \pkg{bamm}.
#' @importFrom methods new
#' @param x An object of class diversity_range

#' @param plot_type Plot type: possible options: "diversity_range"
#' (range-diversity plot),
#'                  "diversity_range_map" (a raster map with
#'                  diversity_range categories),
#'                  "alpha" (a raster map with alpha diversity values),
#'                  "dispersion_field" (a raster with dispersion field)
#' @param legend Logical. If TRUE the legend of the categorical diversity
#' range values will appear.
#' @param legend_position Legend position.
#' @param xlab x label
#' @param ylab y label
#' @param col Plot colors.
#' @param pch Patch type.
#' @param pch_legend Patch type for legends.
#' @param radius Size of the patch for the interactive map.
#' @param ... Graphical parameters. Any argument that can be passed to 1)
#' base::plot, such as axes=FALSE, main='title', ylab='latitude' 2)
#' leaflet::leaflet or 3) leaflet::addCircleMarkers.
#' @return Plot of the results of the diversity_range analysis
#' @details To show interactive diversity_range plots install the 'plotly' R
#' package.
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @rdname plot
#' @export

methods::setMethod(f = "plot",
                   signature = c(x="diversity_range"),
                   function(x,xlab=NULL,plot_type="diversity_range",legend=TRUE,
                            legend_position = "bottomright",
                            ylab=NULL,col=NULL,pch=NULL,pch_legend=19,
                            radius=0.5,...) {
                     if(inherits(x, 'diversity_range')){

                       poptions <- c("diversity_range",
                                     "diversity_range_interactive",
                                     "diversity_range_map",
                                     "alpha","dispersion_field",
                                     "dispersion_field_map")

                       if(!any(poptions %in% plot_type)){
                         stop(paste('plot_type should be one of the following',
                                    'options:',paste0('"',poptions,'"',
                                                     collapse = ", ")))
                       }
                       #slotsin <- methods::slotNames(x)
                       #zeros_alpha <- which(x@alpha ==0)
                       #zeros_disperfield <- which(x@dispersion_field==0)
                       nsites <- x@nsites
                       nsps <- x@nsps

                       if(x@nsites>0){
                         # Random "#000000" = 0
                         # HE/LR "#F6BDC0" = 1
                         # HE/IR  = "#F1A13A" = 2
                         # LE/LR = "#BBDFFA" = 3
                         # HE/HR = #DC1C13" = 4
                         # LE/IR  = "#6987D5" = 6
                         # LE/HR = "#1727AE" = 12
                         codifi <- c("Random" = 0,"HE/LR"=1,"HE/IR"=2,
                                     "LE/LR"=3, "HE/HR"=4, "LE/IR" = 6,
                                     "LE/HR" =12)

                         cols <- c("#000000","#F6BDC0",
                                   "#F1A13A","#BBDFFA",
                                   "#DC1C13","#6987D5",
                                   "#1727AE")

                         names(cols) <- c("Random","HE/LR",
                                          "HE/IR","LE/LR",
                                          "HE/HR","LE/IR",
                                          "LE/HR")

                         COLORES <- cols[names(cols) %in% names(codifi)]
                         #COLORES <- COLORES[c(1:3,5,4,6,7)]
                         if(is.null(pch)){
                           pch <- 19
                         }

                         if(plot_type %in% c("diversity_range",
                              "diversity_range_interactive")){
                           alpha_norm <- x@alpha/nsps
                           dispersion_field <- x@dispersion_field/nsps
                           alpha_st <- x@alpha/x@nsps
                           betty <- round(1/mean(alpha_st),3)
                           fistprom <-x@dispersion_field/x@alpha
                           rho  <-alpha_st*(fistprom-1/betty)
                           am <- min(alpha_st)
                           aM <- max(alpha_st)
                           fm <- min(x@dispersion_field,na.rm = TRUE)/nsites
                           fM <- max(x@dispersion_field,na.rm =TRUE)/nsites

                           rhom  <- min(rho,na.rm = TRUE)
                           rhoM <- max(rho,na.rm = TRUE)

                           vx <- c(am,aM,aM,am)

                           vy <- c(rhom+am/betty,rhom+aM/betty,
                                   rhoM+aM/betty,rhoM+am/betty)

                           if(is.null(col)){
                             col <- x@diversity_range_colors
                           }
                           if(is.null(xlab)){
                             xlab <- expression(alpha)
                           }
                           if(is.null(ylab)){
                             ylab <- "Dispersion field"
                           }
                           xmin1<- 0.97*min(vx)
                           xmax1 <- 1.03*max(vx)
                           ymin1 <- 0.97*min(vy)
                           ymax1 <- 1.03*max(vy)
                           if("diversity_range" == plot_type){
                             plot(alpha_norm,dispersion_field,
                                  xlim=c(xmin1,xmax1),
                                  ylim=c(ymin1,ymax1),
                                  xlab=xlab,ylab=ylab,pch=pch,col=col,...)
                             graphics::lines(graphics::polygon(vx,vy))
                             if(legend){
                               graphics::legend(legend_position,
                                                legend = names(COLORES),
                                                pch =pch_legend,
                                                col = COLORES,bty = "n",...)
                             }
                           }
                           else{
                             Longitude <- x@xy_coordinates[,1]
                             Latitude <-  x@xy_coordinates[,2]
                             labs <- as.factor(x@diversity_range_colors)

                             codifi <- c("Random" = 0,"HE/LR"=1,"HE/IR"=2,
                                         "LE/LR"=3, "HE/HR"=4, "LE/IR" = 6,
                                         "LE/HR" =12)

                             cols <- c("#000000","#F6BDC0",
                                       "#F1A13A","#BBDFFA",
                                       "#DC1C13","#6987D5",
                                       "#1727AE")

                             names(cols) <- c("Random","HE/LR",
                                              "HE/IR","LE/LR",
                                              "HE/HR","LE/IR",
                                              "LE/HR")

                             COLORES <- cols[names(cols) %in% names(codifi)]
                             COLORES <- COLORES[COLORES %in% levels(labs)]

                             levels(labs) <-  names(COLORES)[order(COLORES,
                                                                   levels(labs))
                                                             ]

                             randiv <- x@diversity_range_raster
                             #vals <- stats::na.omit(randiv[])
                             cvals <- raster::cellFromXY(randiv,x@xy_coordinates)
                             vals <- randiv[cvals]
                             cols1 <- ifelse(
                               vals == 0,"#000000",
                               ifelse(vals ==1,
                                      "#F6BDC0",
                                      ifelse(vals==2,
                                             "#F1A13A",
                                             ifelse(vals==3,
                                                    "#BBDFFA",
                                                    ifelse(vals==4,
                                                           "#DC1C13",
                                                           ifelse(vals==6,
                                                                  "#6987D5",
                                            ifelse(vals==12,"#1727AE",NA)))))))

                             div1 <- data.frame(alpha=alpha_norm,
                                                dispersion_field,
                                                Longitude,
                                                Latitude,
                                                labs=as.character(labs),
                                                col=cols1)
                             #div1$col <- cols1

                             #levels(div1$col) <- cols[c(7,2,6,5,4,3,1)]
                             #div1$col <- as.character(div1$col)
                             diversity <- crosstalk::SharedData$new(div1)

                            p1 <-  crosstalk::bscols(
                              plotly::plot_ly() |>
                                  plotly::add_trace(
                                    data=diversity, x = ~alpha,
                                    y = ~dispersion_field,
                                    type="scatter",
                                    mode="markers",
                                    #marker= list(split=cols[c(7,2,6,5,4,3,1)]),
                                    #split = ~labs,
                                    color = ~labs,
                                    colors = COLORES,
                                    inherit = TRUE
                                  )  |>
                                plotly::add_polygons(x=c(vx,vx[1]),
                                                     y=c(vy,vy[1]),
                                                  #name = paste0("Cluster ",1),
                                                     line=list(width=2,
                                                               color="black"),
                                                     fillcolor='transparent',
                                                     hoverinfo = "none",
                                                     showlegend = FALSE,
                                                     inherit = FALSE) |>
                                  plotly::highlight('plotly_selected',
                                                    off = 'plotly_deselect',
                                                    dynamic = FALSE,
                                                    persistent = FALSE),
                              leaflet::leaflet(diversity,height = 900,...) |>
                                leaflet::addTiles() |>
                                leaflet::addCircleMarkers(
                                  radius = radius,
                                  lng = ~Longitude,
                                  lat = ~Latitude,
                                  fillColor = ~col,
                                  color = ~col,opacity = 0.9,...) |>
                                plotly::highlight('plotly_click',
                                                  selectize=FALSE,
                                          off = 'plotly_deselect',
                                          dynamic = FALSE,persistent = FALSE)
                            )
                             print(p1)
                           }

                         }

                         if("alpha" %in% plot_type &&
                            raster::hasValues(x@alpha_raster)){
                           raster::plot(x@alpha_raster,...)
                         }

                         if("dispersion_field" %in% plot_type &&
                            raster::hasValues(x@dispersion_field_raster)){
                           raster::plot(x@dispersion_field_raster,...)
                         }

                         if("diversity_range_map" %in% plot_type &&
                            raster::hasValues(x@diversity_range_raster)){
                           if(is.null(col)){
                             col1 <- cols
                           }else{
                             col1 <- rev(grDevices::terrain.colors(7))
                           }

                           randiv <- x@diversity_range_raster
                           maxval <- raster::maxValue(randiv)
                           r <- raster::ratify(randiv)
                           #rat <- levels(r)[[1]]
                           #rat$VALUE <- c("#000000","#F6BDC0","#F1A13A",
                            #               "#BBDFFA","#DC1C13","#6987D5",
                            #               "#1727AE")
                           #levels(r) <- rat
                           vals <- randiv[]
                           cols1 <- ifelse(
                             vals == 0,"#000000",
                             ifelse(vals ==1,
                                    "#F6BDC0",
                                    ifelse(vals==2,
                                           "#F1A13A",
                                           ifelse(vals==3,
                                                  "#BBDFFA",
                                                  ifelse(vals==4,
                                                         "#DC1C13",
                                                         ifelse(vals==6,
                                                                "#6987D5",
                                                                ifelse(
                                                                  vals==12,
                                                                  "#1727AE",NA)
                                                                ))))))
                           raster::values(r) <- as.factor(cols1)
                           cc <- raster::levels(r)[[1]]
                           raster::plot(r,col=cc$VALUE,legend=FALSE,...)
                           if(legend){
                             graphics::legend(legend_position,
                                              legend = names(COLORES),
                                              pch=15,col = COLORES,
                                              bty = "n",...)
                           }
                         }
                         if("diversity_range_map" %in% plot_type &&
                            !raster::hasValues(x@diversity_range_raster) &&
                            nrow(x@xy_coordinates) == nsites){

                           plot(x@xy_coordinates,col=x@diversity_range_colors,
                                pch=15)
                           graphics::legend("bottomleft",
                                            legend = names(COLORES),
                                            pch=15,col = COLORES,bty = "n")
                         }
                       }

                     }

                   })

#' Predict method of the package \pkg{bamm}.
#' @aliases predict bam-method predict
#' @description predicts species' distribution under suitability changes
#' @param object a of class bam.
#' @param niche_layers A raster or RasterStack with the niche models for
#' each time period
#' @param nbgs_vec A vector with the number of neighbors for the adjacency
#' matrices
#' @param nsteps_vec Number of simulation steps for each time period.
#' @param stochastic_dispersal Logical. If dispersal depends on a probability of
#' visiting neighbor cells (Moore neighborhood).
#' @param disper_prop Probability of dispersal to reachable cells.
#' @param disp_prop2_suitability Logical. If probability of dispersal is
#' proportional
#' to the suitability of reachable cells. The proportional value must be
#' declared
#' in the parameter `disper_prop`.
#' @param animate Logical. If TRUE a dispersal animation on climate change
#' scenarios will be created
#' @param period_names Character vector with the names of periods that will
#' be animated. Default NULL.
#' @param bg_color Color for unsuitable pixels. Default "#F6F2E5".
#' @param suit_color Color for suitable pixels. Default "#0076BE".
#' @param occupied_color Color for occupied pixels. Default "#03C33F".
#' @param ani.width Animation width unit in px
#' @param ani.height Animation height unit in px
#' @param ani.res Animation resolution unit in px
#' @param fmt Animation format. Possible values are GIF and HTML
#' @param filename File name.
#' @param png_keyword A keyword name for the png images generated by
#' the function
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @rdname predict
#' @return A RasterStack of predictions of dispersal dynamics as a function
#' of environmental change scenarios.
#' @examples
#' # rm(list = ls())
#' # Read raster model for Lepus californicus
#' model_path <- system.file("extdata/Lepus_californicus_cont.tif",
#'                           package = "bamm")
#' model <- raster::raster(model_path)
#' # Convert model to sparse
#' sparse_mod <- bamm::model2sparse(model = model,threshold=0.1)
#' # Compute adjacency matrix
#' adj_mod <- bamm::adj_mat(sparse_mod,ngbs=1)
#'
#' # Initial points to start dispersal process
#'
#' occs_lep_cal <- data.frame(longitude = c(-115.10417,
#'                                          -104.90417),
#'                            latitude = c(29.61846,
#'                                         29.81846))
#' # Convert to sparse the initial points
#' occs_sparse <- bamm::occs2sparse(modelsparse = sparse_mod,
#'                                 occs = occs_lep_cal)
#'
#' # Run the bam (sdm) simultation for 100 time steps
#' smd_lep_cal <- bamm::sdm_sim(set_A = sparse_mod,
#'                              set_M = adj_mod,
#'                              initial_points = occs_sparse,
#'                              nsteps = 10)
#' #----------------------------------------------------------------------------
#' # Predict species' distribution under suitability change
#' # scenarios (could be climate chage scenarios).
#' #----------------------------------------------------------------------------
#'
#' # Read suitability layers (two suitability change scenarios)
#' layers_path <- system.file("extdata/suit_change",
#'                            package = "bamm")
#' niche_mods_stack <- raster::stack(list.files(layers_path,
#'                                              pattern = ".tif$",
#'                                              full.names = TRUE))
#' raster::plot(niche_mods_stack)
#' # Predict
#' new_preds <- predict(object = smd_lep_cal,
#'                      niche_layers = niche_mods_stack,
#'                      nsteps_vec = c(50,100))
#'
#' # Generate the dispersal animation for time period 1 and 2
#' \donttest{
#' if(requireNamespace("animation")){
#' ani_prd <- tempfile(pattern = "prediction_",fileext = ".gif")
#' #new_preds <- predict(object = smd_lep_cal,
#' #                      niche_layers = niche_mods_stack,
#' #                      nsteps_vec = c(10,10),
#' #                      animate=TRUE,
#' #                      filename=ani_prd,
#' #                    fmt="GIF")
#'
#'}
#'}

methods::setMethod(f = "predict",
                   signature = methods::signature(object = "bam"),
                   function(object,niche_layers,nbgs_vec=NULL,nsteps_vec,
                            stochastic_dispersal = FALSE,
                            disp_prop2_suitability=TRUE,
                            disper_prop=0.5,
                            animate=FALSE,
                            period_names= NULL,
                            fmt="GIF",filename,
                            bg_color="#F6F2E5",
                            suit_color="#0076BE",
                            occupied_color="#03C33F",
                            png_keyword="sdm_sim",
                            ani.width = 1200,
                            ani.height = 1200,
                            ani.res=300){
                     if(!class(niche_layers) %in% c("RasterLayer",
                                                    "RasterStack",
                                                    "RasterBrick")){
                       stop("niche_layers must be a RasterLayer or RasterStack")
                     }
                     if(methods::is(niche_layers,"RasterLayer"))
                       niche_layers <- raster::stack(niche_layers)
                     n_enm <- raster::nlayers(niche_layers)
                     n_steps_each <- length(nsteps_vec)
                     if(n_steps_each == 1)
                       nsteps_vec <- rep(nsteps_vec,n_enm)
                     if(n_steps_each > 1 && n_enm != n_steps_each){
                       stop(paste("nsteps_vec should have the same",
                                  "length as the number of niche_layers"))
                     }

                     if(!is.null(nbgs_vec)){
                       n_nbgs_each <- length(nbgs_vec)
                       if(n_nbgs_each == 1)
                         nbgs_vec <- rep(nbgs_vec,n_enm)
                       if(n_nbgs_each > 1 && n_enm != n_nbgs_each){
                         stop(paste("nbgs_vec should have the same",
                                    "length as the number of niche_layers"))
                       }
                       ad_mat <- lapply(seq_len(n_nbgs_each), function(x){
                         ad <- bamm::adj_mat(modelsparse = object,
                                             ngbs = nbgs_vec[x])
                       })
                     } else{
                       ad_mat <- lapply(seq_len(n_enm), function(x){
                         methods::new("setM",
                                      adj_matrix=object@adj_matrix,
                                      adj_list = object@adj_list
                                      )
                       })

                     }

                     sim_results <- list(object)

                     initial_points <- object@sdm_sim[[
                       length(object@sdm_sim)]]
                     nsteps <- nsteps_vec[1]
                     niche_mod <- niche_layers[[1]]
                     suit_th <- ifelse(!is.na(object@suit_threshold),
                                       object@suit_threshold,0.1)
                     sparse_mod <- bamm::model2sparse(niche_mod,
                                                      threshold =suit_th)
                     sdm <- bamm::sdm_sim(set_A = sparse_mod,
                                         set_M = ad_mat[[1]],
                                         initial_points = initial_points,
                                         stochastic_dispersal =
                                           stochastic_dispersal,
                                         disp_prop2_suitability =
                                           disp_prop2_suitability,
                                         progress_bar = TRUE,
                                         nsteps = nsteps,
                                         rcpp = TRUE)


                     periods <- raster::stack(object@niche_model,niche_layers)

                     sim_results[[2]] <-sdm
                     if(n_enm > 1){
                       for(x in 2:length(nsteps_vec)){
                         nsteps <- nsteps_vec[x]
                         niche_mod <- niche_layers[[x]]
                         sparse_mod <- bamm::model2sparse(
                           model = niche_mod,
                           threshold = suit_th)
                         bam_object <- sim_results[[x]]
                         initial_points <-
                           bam_object@sdm_sim[[bam_object@sim_steps]]

                         sdm <- bamm::sdm_sim(set_A = sparse_mod,
                                             set_M = ad_mat[[x]],
                                             initial_points = initial_points,
                                             disp_prop2_suitability =
                                               disp_prop2_suitability,
                                             stochastic_dispersal =
                                               stochastic_dispersal,
                                             nsteps = nsteps)
                         sim_results[[x+1]] <- sdm
                       }
                     }

                     names(sim_results) <- paste0("time_period_",
                                                  seq_along(sim_results))

                     if(animate){
                       oldpar <- graphics::par(no.readonly = TRUE)
                       on.exit(graphics::par(oldpar),add=TRUE)
                       nsteps_vec <- c(object@sim_steps,nsteps_vec)
                       nsteps <- sum(nsteps_vec)

                       if(nsteps> 80 && fmt=="GIF"){

                         #contri <- nsteps/length(sim_results)
                         which_steps <- round(seq(1,nsteps,
                                                  along.with = seq_len(80)))
                         step <- max(diff(which_steps))

                        which_stepsL <-  lapply(seq_along(nsteps_vec),
                                                function(x){
                          unique(c(seq(1,nsteps_vec[x],step),
                                   nsteps_vec[x]))
                        })

                       }
                       else{
                         which_stepsL <-  lapply(seq_along(nsteps_vec),
                                                 function(x){
                                                   seq(1,nsteps_vec[x])
                                                 })
                       }

                       titles <- lapply(seq_along(which_stepsL),
                                        function(x) {
                                          if(!is.null(period_names) &&
                                             length(period_names) ==
                                             length(which_stepsL)){
                                            paste0(period_names[x],
                                                   paste0(" (t_"),
                                                   which_stepsL[[x]],")")
                                          }
                                          else{
                                            paste0(paste0("Period ",x," (t_"),
                                                   which_stepsL[[x]],")")
                                          }

                                        })
                       titles <- unlist(titles)


                       sdm_st <- seq_along(sim_results) |>
                         purrr::map(function(x){
                           #nsims <-length(sim_results[[x]]@sdm_sim) -1
                           sdm_ani <- bamm::sim2Raster(
                             sim_results[[x]],
                             which_steps = which_stepsL[[x]])
                           sdm_ani <- sim_results[[x]]@niche_model*1*(sdm_ani+1)

                           return(sdm_ani)
                         })
                       sdm_st <- raster::stack(sdm_st)
                       names(sdm_st) <- titles

                       fmt <- toupper(fmt)
                       if(!fmt %in% c("GIF",'HTML'))
                         stop("fmt should be GIF or HTML")

                       dir1 <- unlist(strsplit(filename,split = "[/]|[\\\\]"))
                       filename <- paste0(dir1,collapse = "/")
                       dir2 <- paste0(dir1[seq_len(length(dir1)-1)],
                                      collapse = '/')
                       if(dir2 == "") dir2 <- "."
                       dir2 <- normalizePath(dir2)
                       if(fmt == "GIF"){
                         animation::ani.options(ani.width = ani.width,
                                                ani.height = ani.height,
                                                ani.res = ani.res)

                         animation::saveGIF({
                           for (i in seq_len(raster::nlayers(sdm_st))) {
                             maxv <- raster::maxValue(sdm_st[[i]])
                             if(maxv<1.5) colores <- c(bg_color,suit_color)
                             else colores <- c(bg_color,suit_color,
                                               occupied_color)

                             graphics::par(xpd = FALSE)


                             raster::plot(sdm_st[[i]],main=titles[i],
                                          col=colores,legend=FALSE,
                                          xaxt = 'n',
                                          yaxt = 'n')

                             graphics::par(xpd = TRUE)
                             graphics::legend(
                               "bottom",
                               legend = c("Unsuitable", "Suitable", "Occupied"),
                               fill = colores,
                               horiz = TRUE,
                               inset = -0.2,
                               cex = 0.75,
                               bty="n"
                             )

                           }

                         },interval=0.8,ani.width = ani.width,
                         movie.name = filename)
                        }
                       if(fmt == "HTML"){

                         dir3 <- file.path(dir2,paste0("pngs_",png_keyword),
                                           fsep = '/')
                         dir3 <- gsub("[\\]","/",dir3)

                         #dir3 <- gsub("[.]","_",dir3)
                         #if(!dir.exists(dir3)) dir.create(dir3)

                         animation::saveHTML({
                           for (i in seq_len(raster::nlayers(sdm_st))) {
                             maxv <- raster::maxValue(sdm_st[[i]])
                             if(maxv<1.5) colores <- c(bg_color,suit_color)
                             else colores <- c(bg_color,suit_color,
                                               occupied_color)

                             graphics::par(xpd = FALSE)


                             raster::plot(sdm_st[[i]],main=titles[i],
                                          col=colores,legend=FALSE,
                                          xaxt = 'n',
                                          yaxt = 'n')

                             graphics::par(xpd = TRUE)
                             graphics::legend(
                               "bottom",
                               legend = c("Unsuitable", "Suitable", "Occupied"),
                               fill = colores,
                               horiz = TRUE,
                               inset = -0.2,
                               cex = 0.75,
                               bty="n"
                             )
                           }
                         },img.name = png_keyword,
                         imgdir = dir3 ,
                         htmlfile = filename,
                         ani.width=ani.width,
                         ani.height=ani.width,interval=0.1,
                         ani.dev = function(...){
                           grDevices::png(res=ani.res,...)})
                       }


                     }

                     return(sim_results)
                   })



