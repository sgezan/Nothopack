#' Predicts and/or projects stand basal area for other species
#'
#' Evaluates stand-level input variables to predict or project basal area for other species.
#' Projections are based in 1 year increments.
#'
#' @param AD0 Dominant age (years) of the stand at current time
#' @param BA990 Basal area (m2/ha) of other species at current time (required for projections)
#' @param PNHAN0 Proportion of number of trees of Nothofagus (trees/ha) of the stand at time 0
#' @param PNHAN1 Proportion of number of trees of Nothofagus (trees/ha) of the stand at time 1
#' @param PBAN0 Proportion of basal area (m2/ha) of Nothofagus of the stand at time 0
#' @param PBAN1 Proportion of basal area (m2/ha) of Nothofagus of the stand at time 1
#' @param projection if TRUE projection from BA0 is executed for a 1 year increment
#'
#' @return Basal area for other species (BA99, m2/ha) for the current age (for prediction)
#' or at AD0+1 (for projection)
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example 1: Predicts Basal Area
#' BA99est <- BA99module(AD0=19, PNHAN0=0.81, PBAN0=0.91, projection=FALSE)
#' BA99est$BA990
#'
#' # Example 2: Projects Basal Area
#' BA99est <- BA99module(BA990=3.2, AD0=19, PNHAN0=0.81, PNHAN1=0.81,
#'                       PBAN0=0.91, PBAN1=0.91, projection=TRUE)$BA991
#' BA99est

BA99module <- function(BA990=NA, AD0=NA, PNHAN0=NA, PNHAN1=NA, PBAN0=NA, PBAN1=NA, projection=FALSE){

  # Model 1 (linear): BA99 = exp(b0)*EDOM^b1*PNHAN^b2*PBAN^b3 + c
  bm<-c(2.04068, 0.07997, -0.13080, -2.22934, -10) # b0,b1,b2,b3,c  - Old set of parameters substracts 10
  bm<-c(1.99503, 0.09436,	-0.21578, -1.87264, 0) # b0,b1,b2,b3,c  - Paremeters form manuscript

  # Prediction
  if (!projection){
    BA990<-exp(bm[1])*(AD0^bm[2])*(PNHAN0^bm[3])*(PBAN0^bm[4]) + bm[5]
    BA991<-NA
  }

  # Projection
  if (projection){
    BA991<-BA990*exp(bm[2]*log((AD0+1)/AD0) + bm[3]*log(PNHAN1/PNHAN0) + bm[4]*log(PBAN1/PBAN0))
  }

  return(list(BA990=BA990, BA991=BA991))
}
