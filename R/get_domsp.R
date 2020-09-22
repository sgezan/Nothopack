#' Determines the dominant specie of a Nothofagus stand
#'
#' It determines the dominant Nothofagus specie of a given stand
#' based on stand-level basal area.
#'
#' @param BA Vector of basal area (m2/ha) of the stand (1: Rauli, 2: Roble, 3: Coigue, 4:Others)
#'
#' @return DOM.SP Dominant specie of the stand (1: Rauli, 2: Roble, 3: Coigue, 4:Mixed)
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065, Chile
#'
#' @examples
#' # Example: Obtaining species from vector of Basal Area
#' BA <- c(4.5,2.8,41.6,2.4)
#' (DOM.SP <- get_domsp(BA))

get_domsp <- function(BA=NA){

  # Proportion of SPECIES by basal area
  BA0 <- sum(BA, na.rm = TRUE)
  BAN <- sum(BA[1:3], na.rm = TRUE)   # BA for all Nothofagus
  PBAN <- BAN/BA0   # Proportion BA for all Nothofagus
  PBA1 <- BA[1]/BAN   # Rauli
  PBA2 <- BA[2]/BAN   # Roble
  PBA3 <- BA[3]/BAN   # Coigue
  PBA99 <- BA[4]/BA0   # Others
  Proportion<-c(PBA1,PBA2,PBA3)

  # Obtaining dominant SPECIES.
  if (PBAN < 0.6 ){    # If BA Nothofagus represent less than 60% of the stand
    warning('Stand is not dominated by Nothofagus (PBAN < 0.6). Returning plot with highest proportion of Basal Area')
    POS.MAX  = which(order(Proportion,decreasing = TRUE)==1)[1]  #Position of the bigest condition
    DOM.SP <- POS.MAX

    # DOM.SP <- 99  # Others besides Nothofagus
  } else {
    if (PBA1 >= 0.7){
      DOM.SP <- 1  # Rauli
    } else if (PBA2 >= 0.7){
      DOM.SP <- 2  # Roble
    } else if (PBA3 >= 0.7){
      DOM.SP <- 3  # Coigue
    } else {
      #If it is mixed, return the Nothofagus species with highest PBA
      POS.MAX  = which(order(Proportion,decreasing = TRUE)==1)[1]  #Position of the bigest condition
      DOM.SP <- POS.MAX
      warning('There is no individual Nothofagus specie with a proportion of Basal Area > 70%')

      }
  }

  return(DOM.SP)
}


