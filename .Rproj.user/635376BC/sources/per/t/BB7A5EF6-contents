# R classes for bamm package
# December 2019
# Version 0.3.2
# Licence GPL v3

#' S4 classes to organize data and results of \code{bamm} objects
#' @aliases coordinates-class
#' @importFrom methods new
#' @slot coordinates A two column matrix with coordinates
#' @slot eigen_vec Eigen vector of adjacency matrix
#' @slot eigen_val Eigen value of adjacency matrix
#' slot g_model A raster representing the geographic area
#' slot g_sparse A sparse matrix of the geographic area
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @return An object of class g_area
#' @exportClass g_area
#' @export
#'
g_area <- methods::setClass(Class = "g_area",
                            slots = c(coordinates="matrix",
                                      eigen_vec = "matrix",
                                      eigen_val = "numeric"
                                      #g_raster = "RasterLayer",
                                      #g_sparse = "dgCMatrix"
                                      ))

#' Class for the A set of the BAM diagram
#'
#' A class for the A set of the BAM diagram. It contains raster models and
#' IDs of pixels with values different than NA.
#'
#' @slot niche_model A niche model in raster format. It can be a binary model
#' or continuous.
#' If the model is in a continuous format.
#' @slot suit_threshold Suitability value used to binarize continuous model
#' @slot cellIDs A numeric vector with the IDs of the cells with prediction
#' values
#' @slot suit_values A numeric vector with suitability value of the continuous
#' map
#' @slot sparse_model A niche model in sparse matrix format
#' @importClassesFrom raster RasterLayer
#' @importClassesFrom raster RasterStack
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods new
#' @aliases A-class
#' @return An object of class setA
#' showClass("setA")
#' @rdname setA
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @export
#' @exportClass setA
#'
setA <- methods::setClass(Class = "setA",
                          slots = c(niche_model = "RasterLayer",
                                    suit_threshold = "numeric",
                                    cellIDs = "numeric",
                                    suit_values = "numeric",
                                    sparse_model = "dgCMatrix"),
                          contains = "g_area",
                          validity = function(object){
                            model_class <- class(object@niche_model)
                            if(!model_class %in% c("RasterLayer",
                                                 "RasterStack")){
                              stop("niche_model should be a raster layer
                                   or raster stack")
                            }
                            })

#' Class for the M set of the \code{bamm} diagram
#' @aliases M-class
#' @slot adj_matix An adjacency matrix
#' @slot adj_list An adjacency list
#' @slot initial_points A presence-absence vector with species' occurrences
#' @slot n_initial_points Number of initial points used to start the dispersal
#' process
#' @slot ngbs Number of neighbors
#' @importClassesFrom raster RasterLayer
#' @importClassesFrom raster RasterStack
#' @importClassesFrom Matrix dgCMatrix
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @return An object of class setM
#' @examples
#' showClass("setM")
#' @exportClass setM

setM <- methods::setClass(Class = "setM",
                          contains = "g_area",
                          slots =  c(adj_matrix = "dgCMatrix",
                                     adj_list = "list",
                                     initial_points = "dgCMatrix",
                                     ngbs = "numeric"))

#' Class \code{bam} digram
#' @aliases bam-class
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot sdm_sim A list of sparse vectors representing the area occupied
#' @slot palatable_matrices A list of sparse vectors representing palatable
#' sites.
#' @slot sim_steps Number of simulation steps
#' by the species
#' @return An object of class bam
#' @exportClass bam

bam <- methods::setClass(Class = "bam",
                         contains = c("setA","setM"),
                         slots = c(sdm_sim = "list",
                                   palatable_matrices = "list",
                                   sim_steps="numeric"))


#' Class  \code{community_sim} digram
#' @aliases community-class
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot community_sim A list of sparse vectors representing the area occupied
#' by the species
#' @return An object of class community_sim
#' @exportClass community_sim

community_bam <- methods::setClass(Class = "community_sim",
                                   #contains = c("bam"),
                                   slots = c(community_sim = "list"))


#' Class  \code{pam} Presence-Absence Matrix
#' @aliases PAM
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot pams A list of sparse matrices representing Presence-Absence Matrix for
#' each simulation time
#' @slot which_steps Simulation steps
#' @slot sp_names Names of species in the PAM
#' @slot grid Raster grid of the studied area
#' @slot cellIDs Cells with ids of the PAM sites
#' @exportClass pam
#' @return An object of class pam
#' @export
#'
pam <- methods::setClass(Class = "pam",
                         #contains = c("bam"),
                         slots = c(pams = "list",
                                   which_steps= "numeric",
                                   sp_names = "character",
                                   grid = "RasterLayer",
                                   cellIDs = "numeric"))



#' Class  \code{leaflet} leaflet
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @return An object of class leaflet
leaflet <- methods::setClass(Class = "leaflet")

#' Class  \code{csd}
#' @aliases csd-class
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot connections A data.frame with four columns: x, y, clusterID and
#' cluster_size
#' @slot interactive_map A leaflet map with markers showing the geographical
#' clusters
#' @slot raster_map A raster map with cluster IDs as values.
#' @exportClass csd
#' @export
#' @return An object of class csd
#' @importFrom leaflet leaflet
csd <- methods::setClass(Class = "csd",
                         representation(connections="data.frame",
                                        interactive_map = "leaflet",
                                        raster_map = "RasterLayer"))




#setClassUnion("sparse_or_matrix", c("dgeMatrix", "matrix"))


#' Class  \code{bioindex_sparse}
#' @aliases biodiversity_index_sparse
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot alpha A sparse matrix with the richness of species per site
#' @slot omega A sparse matrix with the range size of every species
#' @slot wBeta A numeric value with Whittaker’s multiplicative beta index
#' @slot laBeta A numeric value with Lande’s additive beta index
#' @slot leBeta A numeric value with Legendre’s beta index
#' @slot nestedness A numeric value with Wright and Reeves' nestedness
#' @slot dispersion_field A sparse matrix with the set of ranges of all species
#' that occur in at each locality
#' @slot richness_field A sparse matrix with the number of  shared species in
#' each site
#' @exportClass bioindex_sparse
#' @return An object of class bioindex_sparse
#' @export

bioindex_sparse <- methods::setClass(Class = "bioindex_sparse",
                              representation (
                                alpha = "dgeMatrix",
                                omega = "dgeMatrix",
                                wBeta = "numeric",
                                laBeta = "numeric",
                                leBeta = "numeric",
                                nestedness = "numeric",
                                dispersion_field = "dgeMatrix",
                                richness_field = "dgeMatrix"
                              ))


#' Class  \code{bioindex}
#' @aliases biodiversity_index
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot alpha A matrix with the richness of species per site
#' @slot omega A matrix with the range size of every species
#' @slot wBeta A numeric value with Whittaker’s multiplicative beta index
#' @slot laBeta A numeric value with Lande’s additive beta index
#' @slot leBeta A numeric value with Legendre’s beta index
#' @slot nestedness A numeric value with Wright and Reeves' nestedness
#' @slot dispersion_field A matrix with the set of ranges of all species that
#' occur in at each locality
#' @slot richness_field A matrix with the number of  shared species in
#' each site
#' @exportClass bioindex
#' @return An object of class bioindex
#' @export

bioindex <- methods::setClass(Class = "bioindex",
                                     representation (
                                       alpha = "matrix",
                                       omega = "matrix",
                                       wBeta = "numeric",
                                       laBeta = "numeric",
                                       leBeta = "numeric",
                                       nestedness = "numeric",
                                       dispersion_field = "matrix",
                                       richness_field = "matrix"
                                     ))


#' Class  \code{diversity_range}
#' @aliases diversityrange
#' @author Luis Osorio-Olvera & Jorge Soberón
#' @slot alpha A column vector with species richness per site
#' @slot omega A column vector with the size of the area of distribution per
#' species.
#' @slot alpha_raster Species richness in raster format.
#' @slot dispersion_field A matrix with the set of ranges of all species that
#' occur in at each locality.
#' @slot dispersion_field_raster Raster object with the observed values of
#' dispersion field.
#' @slot diversity_range_raster Raster object of diversity range.
#' @slot diversity_range_colors Colors to plot endemism levels.
#' @slot null_dispersion_field_dist A matrix with dispersion field null
#' distribution.
#' @slot xy_coordinates A matrix of geographical coordinates
#' @slot n_iterations Number of iterations used to estimate the dispersion
#' field null distribution.
#' @slot nsps Number of species in the PAM.
#' @slot nsites Number of sites in the PAM.
#' @return An object of class diversity_range
#' @exportClass diversity_range
#' @export

diversity_range <- methods::setClass(Class = "diversity_range",
                              representation (
                                alpha = "matrix",
                                omega = "matrix",
                                alpha_raster = "RasterLayer",
                                dispersion_field = "matrix",
                                dispersion_field_raster ="RasterLayer",
                                diversity_range_raster ="RasterLayer",
                                diversity_range_colors ="character",
                                null_dispersion_field_dist ="matrix",
                                xy_coordinates = "matrix",
                                n_iterations = "numeric",
                                nsps ="numeric",
                                nsites = "numeric"
                              ))
