#' range_diversity_analysis: diversity analysis
#' @description diversity_range_analysis biodiversity indices related to
#' diversity-range plots
#' @param pam A Presence-Absence-Matrix of matrix class or sparse matrix.
#' @param xy_mat A two dimensional matrix with longitude and latitude data.
#' @param lower_interval Lower interval.
#' @param upper_interval Upper interval.
#' @param raster_templete A raster template.
#' @param niter Number of iterations to obtain the distribution.
#' @param return_null_dfield If TRUE the null distribution of dispersal field
#' will be returned.
#' @param randal Randomization algorithm applied to the PAM.
#' Possible choices "curveball","fastball" and "indep_swap".
#' @param parallel If TRUE the computations will be performed in parallel.
#' @param n_cores Number of cores for the parallel computation.
#' @return An object of class  \code{\link[bamm]{diversity_range}}. The main
#' result is the diversity range analysis which shows jointly two indices
#' describing the community composition of every cell in the grid:
#' (1) the relative number of species, and (2) the mean dispersion field
#' (see plot method for \code{\link[bamm]{plot}} (Soberon et al. 2022).
#' The contains 12 slots with different measurements of biodiversity such
#' as alpha diversity (species richness in each site or pixel),
#' omega (size of the area of distribution of each species),
#' dispersion field (the standardized size of the area of distribution of
#' all species occurring in each pixel).
#' @details For more information about the biodiversity indices
#' see Soberon and Cavner (2015). For detail about the diversity range analysis
#' see Soberon et al. (2022). To plot diversity range results use
#' \code{\link[bamm]{plot}} method for objects of class
#' \code{\link[bamm]{diversity_range}}.
#'
#' For details about randomization algorithms applied to the PAM see
#' \code{\link[bamm]{null_dispersion_field_distribution}}.
#'
#' @references
#' \insertRef{Soberon2022}{bamm}.
#'
#' \insertRef{Soberon2015}{bamm}.
#'
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @importFrom future plan tweak multisession sequential
#' @importFrom Rcpp sourceCpp
#' @useDynLib bamm
#' @export
#' @examples
#' set.seed(111)
#' pam <- matrix(rbinom(10000,1,0.5),nrow = 100,ncol = 1000)
#' rdivan <- bamm::diversity_range_analysis(pam=pam,
#'                                          parallel = FALSE,
#'                                          niter = 10,
#'                                          return_null_dfield=TRUE)
#' bamm::plot(rdivan,plot_type="diversity_range")
#' # Lagomorphos
#' \donttest{
#' lagos_path <- system.file("extdata/conejos",
#'                           package = "bamm")
#' enm_path <- list.files(lagos_path,
#'                        pattern = ".tif",
#'                        full.names = TRUE)

#' en_models <- raster::stack(enm_path) >0.01
#' nonas <- which(!is.na(en_models[[1]][]))
#' xy_mat <- sp::coordinates(en_models[[1]])[ nonas,]
#' pam <- bamm::models2pam(en_models,sparse=FALSE)
#'
#' rdivan <- bamm::diversity_range_analysis(pam=pam,
#'                                          xy_mat=xy_mat,
#'                                          raster_templete = en_models[[1]],
#'                                          parallel=TRUE,
#'                                          n_cores=2,
#'                                          return_null_dfield=TRUE)
#' bamm::plot(rdivan,plot_type="diversity_range")
#' bamm::plot(rdivan,plot_type="diversity_range_map")
#' if(require(plotly) && require(crosstalk)){
#' #bamm::plot(rdivan,plot_type="diversity_range_interactive")
#' }
#' }

diversity_range_analysis <- function(pam,xy_mat=NULL,lower_interval=0.05,
                                     upper_interval=0.95, raster_templete=NULL,
                                     niter=100,return_null_dfield=FALSE,
                                     randal = "indep_swap",
                                     parallel=TRUE,
                                     n_cores=2){

  ral <- match.arg(arg = randal,
                   choices = c("indep_swap","curveball","fastball"))

  if(!methods::is(pam,"matrix") & !is.numeric(pam[1,1])){
    stop("pam object should be a binary matrix")
  }

  results <- methods::new(Class = "diversity_range")

  distfield_rand <- null_dispersion_field_distribution(pam = pam,
                                                       n_iter=niter,
                                                       parallel=parallel,
                                                       randal = randal,
                                                       n_cores =n_cores)

  if(return_null_dfield){
    results@null_dispersion_field_dist <- distfield_rand
  }


  bioind <- bamm::pam2bioindex(pam=pam,
                               biodiv_index = c("alpha",
                                                "dispersion_field"),
                               as_sparse = FALSE)

  results@alpha <- bioind@alpha
  results@omega <- bioind@omega
  results@nsps <- ncol(pam)
  results@nsites <- nrow(pam)
  results@n_iterations <- niter
  results@dispersion_field <- bioind@dispersion_field
  disfield_cat <-  null_dispersion_field_cat(dfield = bioind@dispersion_field,
                                             dfield_rand = distfield_rand,
                                             lower_interval=lower_interval,
                                             upper_interval=upper_interval)


  #cols <- c("#000000","#F6BDC0",
  #          "#F07470","#BBDFFA",
  #         "#DC1C13","#6987D5",
  #          "#1727AE")
  cols <- c("#000000","#F6BDC0",
            "#F1A13A","#BBDFFA",
            "#DC1C13","#6987D5",
            "#1727AE")
  names(cols) <- c("Random","HE/LR",
                   "HE/IR","LE/LR",
                   "HE/HR","LE/IR",
                   "LE/HR")

  richness_cat <-  bioind@alpha
  q_rich <-  stats::quantile(richness_cat)
  q_rich_low <- q_rich[2]
  q_rich_int <- q_rich[4]
  q_rich_hig <- q_rich[5]

  q_rich_low_ids <- which(richness_cat[,1] <  q_rich_low)
  q_rich_int_ids <- which(richness_cat[,1] >= q_rich_low &
                            richness_cat[,1] < q_rich_int)
  q_rich_hig_ids <- which(richness_cat[,1] >= q_rich_int &
                            richness_cat[,1] <= q_rich_hig)

  richness_cat[q_rich_low_ids,1] <- 1
  richness_cat[q_rich_int_ids,1] <- 2
  richness_cat[q_rich_hig_ids,1] <- 4
  dfalpha <- disfield_cat*richness_cat
  range_div_cols <-  dfalpha
  #alpha_rasterC <-   range_div_cols
  range_div_cols <- as.factor(dfalpha)

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
  levels(range_div_cols) <- codifi[codifi %in% levels(range_div_cols)]
  vals <- as.numeric(as.character(range_div_cols))
  #levels(range_div_cols) <- cols[codifi %in% levels(range_div_cols)]
  cl <- function(vals){
    ifelse(vals == 0,"#000000",
           ifelse(vals ==1, "#F6BDC0",
                  ifelse(vals==2,"#F1A13A",
                         ifelse(vals==3,"#BBDFFA",
                                ifelse(vals==4,"#DC1C13",
                                       ifelse(vals==6,"#6987D5",
                                              ifelse(vals==12,"#1727AE",
                                                     NA)))))))
  }
  results@diversity_range_colors <- cl(vals = vals)
  if(is.matrix(xy_mat) || is.data.frame(xy_mat)){


    results@xy_coordinates <- data.matrix(xy_mat)

    if(methods::is(raster_templete,"RasterLayer")){

      r1 <- raster_templete*0
      alpha_raster <- r1
      names(alpha_raster) <- "alpha"
      dispersion_field_raster <- r1
      names(dispersion_field_raster) <- "dispersion_field"
      diversity_range_raster <- r1
      names(diversity_range_raster) <- "diversity_range"
      cellIDs <- raster::cellFromXY(r1,xy_mat[,1:2])
      alpha_raster[cellIDs]<- bioind@alpha
      results@alpha_raster <- alpha_raster
      dispersion_field_raster[cellIDs] <- bioind@dispersion_field
      diversity_range_raster[cellIDs] <- vals
      results@dispersion_field_raster <- dispersion_field_raster
      results@diversity_range_raster<- diversity_range_raster
    }

  }

  return(results)

}
