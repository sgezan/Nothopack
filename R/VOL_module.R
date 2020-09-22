#' Estimates stand-level volume based on stand-level parameters
#'
#' It reads stand-level input variables (BA, HD and PNHAN, PBAN) to estimate stand-level total
#' volume with bark (m3.ha). There are two models: 1:requires BA, HD, PNHAN and PBAN, and
#' 2:requires BA, HD, PNHAN. The later model is more stable.
#'
#' @param BA Basal Area (m2/ha) of the current stand
#' @param HD Dominant height (m) of the current stand
#' @param PNHAN Proportion of Nothofagis trees in the stand (optional)
#' @param PBAN Proportion of Nothofagus basal area in the stand (optional)
#'
#' @return Total volume with bark (m3/ha) for the current conditions across all species
#'
#' @author
#' S.A. Gezan, P. Moreno, S. Palmas
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#' (updated version with recalculated stand volumes)
#'
#' @examples
#' # Example 1: Predicts volume with PNHAN and PBAN
#' Vmodule(BA=32, HD=15.9, PNHAN=0.81, PBAN=0.90)
#'
#' # Example 2: Predicts volume with PNHAN
#' Vmodule(BA=32, HD=15.9, PNHAN=0.81)
#'
#' # Example 3: Predicts volume
#' Vmodule(BA=32, HD=15.9)

Vmodule <- function(BA=NA, HD=NA, PNHAN=NA, PBAN=NA){

  # Any BA or HD are missing
  if(is.na(BA) || is.na(HD)){
    stop('Missing BA or HD for volume estimation in VOL_module')
  }

  # Updated version with recalculated stand volumes and Baskerville correction
  if (!is.na(PNHAN)) {
    if (!is.na(PBAN)) {
      VOL <- exp(-0.934+0.03576/2)*(BA^0.9100)*(HD^1.0949)*(PNHAN^-0.270)*(PBAN^0.2241)
    } else {
      VOL <- exp(-0.820+0.03673/2)*(BA^0.9076)*(HD^1.0602)*(PNHAN^0.1197)
    }
  } else {
    VOL <- exp(-0.735+0.03891/2)*(BA^0.8662)*(HD^1.0661)
  }

  return(VOL)
}
