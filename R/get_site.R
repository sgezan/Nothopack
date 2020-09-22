#' Calculates missing stand-level variable from the set: AD, HD, SI
#'
#' When two of stand level variables: dominat age (AD),
#' dominant height (HD) and site index (SI) are given, it returns the value of
#' the remaining stand-level parameter.
#'
#' @param DOM.SP Dominant species (1: Rauli, 2: Roble, 3: Coigue)
#' @param ZONE Growth zone (1, 2, 3, 4)
#' @param AD Dominant age (year) of the stand
#' @param HD Dominant height (m) of dominant specie in the current stand
#' @param SI Site index (m) at reference dominant age of 20 years of the stand
#'
#' @return The missing stand-level parameter: AD, HD or SI
#'
#' @author
#' S. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065, Chile
#'
#' Gezan, S.A. and Moreno, P. (2000). Curvas de Sitio - Altura dominante para renovales
#' de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065, Chile
#'
#' @examples
#' # Example 1: Obtain Dominant Age
#' (AD <- get_site(DOM.SP=1, ZONE=2, HD=14, SI=10))
#'
#' # Example 2: Obtain Dominant Height
#' (HD <- get_site(DOM.SP=1, ZONE=2, AD=25, SI=10))
#'
#' # Example 3: Obtain Site Index
#' (SI <- get_site(DOM.SP=1, ZONE=2, AD=19, HD=13.5))

get_site <- function(DOM.SP, ZONE, AD=NA, HD=NA, SI=NA){

  # Correct Model is: HD = a [1 – {1 – (IS / a) c } ((E - 2) / (20 - 2)] 1/c
  #                   c = b0 + b1 IS

  coef.list <- subset(hd_coef, hd_coef_zone == ZONE & hd_coef_sp_code == DOM.SP,
                      select = c(hd_coef_a, hd_coef_b0, hd_coef_b1) )

  if (sum(is.na(c(AD, HD, SI))) >= 2 ){
    parm<-NA
    warning('There must be at least two stand parameters provided')

  } else if ( sum(is.na(c(AD, HD, SI))) == 0 ){
    warning('Why would you use this calculation? You have already have the three variables')

  } else if (is.na(AD)) {               # If Initial age is missing

    # Definition of AD function for bisection method (should it be somewhere else?)
    AD.eq <- function(x){
      c <- coef.list$hd_coef_b0 + coef.list$hd_coef_b1*SI
      - HD + coef.list$hd_coef_a * (1-(1-(SI/coef.list$hd_coef_a)^c)^((x-2)/18) )^(1/c)
    }
    # Bisection method
    #parm <- tryCatch(pracma::bisect(AD.eq,
    #                                a=2, b=100,
    #                                maxiter=100)$root)

    # Brute search method
    seqX <- seq(1,90,0.01)
    lista1 <- unlist(lapply(X = seqX, FUN = AD.eq))   #lista de resultados.
    posit <- which(abs(lista1 - 0) == min(abs(lista1 - 0), na.rm = TRUE))
    parm <- seqX[posit]


  } else if (is.na(HD)) {              # If Dominant height is missing
    c <- coef.list$hd_coef_b0 + coef.list$hd_coef_b1*SI
    parm <- coef.list$hd_coef_a * (1-(1-(SI/coef.list$hd_coef_a)^c)^((AD-2)/18))^(1/c)
    #return(parm)

  } else if (is.na(SI)) {              # If Site Index is missing

    # Definition of SI function for bisection method (should it be somewhere else?)
    # Since the function requires the coefficients and I think the bisect method
    # looks for x in the function, it is probably safer to define it here.
    SI.eq <- function(x){
      c <- coef.list$hd_coef_b0 + coef.list$hd_coef_b1 * x
      - HD + coef.list$hd_coef_a * (1-(1-(x/coef.list$hd_coef_a)^c)^((AD-2)/18))^(1/c)
    }

    # Brute search method
    seqX <- seq(1,50,0.01)
    lista1 <- unlist(lapply(X = seqX, FUN = SI.eq))   #lista de resultados. Esto solo busca IS de 1 a 40
    posit <- which(abs(lista1 - 0) == min(abs(lista1 - 0), na.rm = TRUE))
    parm <- seqX[posit]

  }

  return(round(parm,2))

}

# Note
#   - Need to restict input of SI (0-40), AD (2-100) to reasonable numbers

