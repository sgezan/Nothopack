#' Calculates missing stand-level variable from the set: N, BA, QD
#'
#' When two of the stand level variables: basal area (BA, m2/ha),
#' number of trees (N, trees/ha) or quadratic diameter (QD, cm) are given, it
#' returns the value of the remaining missing stand-level parameter.
#'
#' @param BA Basal area (m2/ha) of the stand
#' @param N Number of trees (trees/ha) of the stand
#' @param QD Quadratic diameter (cm) of the stand
#'
#' @return The missing stand level parameter: BA, N or QD
#'
#' @author
#' S. Palmas, S.A. Gezan and P. Moreno
#'
#' @examples
#' # Example 1: Obtain quadratic diameter
#' (QD <- get_stand(BA=25.2, N=1400))
#'
#' # Example 2: Obtain basal area
#' (BA <- get_stand(QD=17, N=1400))
#'
#' # Example 3: Obtain number of trees
#' (N <- get_stand(BA=25.2, QD=17))

get_stand <- function(BA=NA, N=NA, QD=NA){

  if (sum(is.na(c(BA, QD, N))) >= 2 ){
    stop('There must be at least two stand parameters provided')

  } else if ( sum(is.na(c(BA, QD, N))) == 0 ){
    warning('Why would you use this calculation? You already have the three variables')

  } else if (is.na(BA)) {      # Estimation of BA witn N and QD
    BA <- (pi/4)*N*(QD/100)^2
    return(BA)

  } else if (is.na(N)) {       # Estimation of N with BA and QD
    if (QD == 0){
      N <- 0
    } else { N <- (4/pi)*BA*(100/QD)^2 }
    return(N)

  } else if (is.na(QD)) {      # Estimation of QD with N and BA
    if (N == 0 | BA == 0){
      QD <- 0
    } else { QD <- 100*((4/pi)*(BA/N))^0.5 }
    return(QD)
  }
}
