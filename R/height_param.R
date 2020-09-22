#' Estimates total individual tree height in function of its DBH and stand-level parameters
#'
#' Estimates individual tree height form a model of stand-level parameters:
#' dominant height (HD, m), quadratic diameter (QD, cm), and diameter at breast height (DBH, cm).
#'
#' @param DOM.SP Dominant species (1: Rauli, 2: Roble, 3: Coigue) of the stand
#' @param ZONE Growth zone (1, 2, 3, 4) of the stand
#' @param HD Dominant height (m) of dominant specie in the current stand
#' @param QD Quadratic diameter (cm) of the stand
#' @param DBH Diameter at breast height (cm) of tree
#'
#' @return Individual total tree height HT (m)
#'
#' @author
#' S. Palmas, S.A. Gezan and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example: Estimates HT for trees with DBH=24 for different stands
#' (HT <- height_param(DOM.SP=2, ZONE=2, HD=15, QD=12, DBH=24))
#' (HT <- height_param(DOM.SP=2, ZONE=2, HD=15, QD=18, DBH=24))
#' (HT <- height_param(DOM.SP=2, ZONE=2, HD=19, QD=12, DBH=24))
#' (HT <- height_param(DOM.SP=2, ZONE=2, HD=19, QD=18, DBH=24))

height_param <- function(DOM.SP, ZONE, HD=NA, QD=NA, DBH=NA){

  if(is.na(DOM.SP) || is.na(ZONE) || is.na(HD) || is.na(QD) || is.na(DBH)){
    warning('Missing input information for height_param module.')
  }

  coef.list <- subset(hparam_coef, hparam_zone == ZONE & hparam_dom_sp_code == DOM.SP,
                      select = c(hparam_b0, hparam_b1, hparam_b2, hparam_b3, hparam_b4, hparam_b5))

  hest <-(coef.list$hparam_b0 + coef.list$hparam_b1*HD + coef.list$hparam_b2*(QD^0.95)
  + coef.list$hparam_b3*exp(-0.08*DBH) + coef.list$hparam_b4*(HD^3)*exp(-0.08*DBH)
  + coef.list$hparam_b5*(QD^3)*exp(-0.08*DBH))

  return(hest)
}
