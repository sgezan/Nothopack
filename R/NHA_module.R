#' Projects number of trees (N) from stand-level parameters to the next year
#'
#' Evaluates stand-level input variables to project number of trees per hectare
#' to the following year. This is a mortality model, as it does not incorporate ingrowth.
#'
#' @param NHA0 Number of trees (trees/ha) of the stand at current time 0
#' @param QD0 Quadratic diameter (cm) of the stand at current time 0
#' @param DOM.SP Dominant specie (1:Rauli, 2:Roble, 3:Coigue, 4:Others or Mixed)
#' @param model Number of fitted model for N estimation (1:Original Reineke,
#'              2:Re-fitted Reineke, 3:Reineke with correction factor) (default = 3)
#' @param return.DQmax If true returns DQ maximum from Reineke's model (default = FALSE)
#'
#' @return Number of trees (N1, m2/ha) for the next year (i.e. time 1)
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' Gezan, S.A., Ortega, A., Andenmatten, E. (2007). Diagramas de manejo de densidad para renovales
#' de roble, rauli y coigue en chile. Bosque (Valdivia) 28, 97-105
#'
#' @examples
#' # Example: Number of trees for next year with current NHA0
#' (N1 <- NHAmodule(NHA0=2730, QD0=12.43, model=1))
#' (N1 <- NHAmodule(NHA0=2730, QD0=12.43, model=2))

NHAmodule <- function(NHA0=NA, QD0=NA, DOM.SP=4, model=3, return.DQmax=FALSE){

  # Model: log(n_trees_ha2) = log(n_trees_ha1)*(1 - theta*Delta.ANHO*(dq1/dq_max_original))
  if (model == 1){

    intercepts <- c(11.6167, 11.3770, 11.7630, 11.6167) #From GOA2007 Rauli, Roble, Coigue, Mixto
    int.SP <- intercepts[DOM.SP]
    theta <- 0.003595746 # Using original Reineke function and esimated theta (from Palmas manuscript)
    DQmax <- exp((log(NHA0) - int.SP)/-1.4112)
    NHA1 <- exp(log(NHA0)*(1-theta*(QD0/DQmax)))

  } else if (model == 2) {    #Modelo no muy revisado and only with one species

    theta <- 0.0056560  # Using new Reineke fitted function and estimated theta
    DQmax <- exp((log(NHA0) - 13.500416)/-1.990455)
    NHA1 <- exp(log(NHA0)*(1-theta*(QD0/DQmax)))

  } else if (model == 3) {    # Original model with bias correction k - based on remeasured plots

    k <- 1.014
    intercepts <- c(11.6167, 11.3770, 11.7630, 11.6167) #From GOA2007 Rauli, Roble, Coigue, Mixto
    int.SP <- intercepts[DOM.SP]
    theta <- 0.003595746 # Using original Reineke function and esimated theta (from Palmas manuscript)
    DQmax <- exp((log(NHA0) - int.SP)/-1.4112)
    NHA1 <- k*exp(log(NHA0)*(1-theta*(QD0/DQmax)))

  }

  if(return.DQmax){
    return(DQmax)
  } else {
    return(NHA1)
  }
}
