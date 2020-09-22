#' Calculates stem diameter or height based on taper equation models
#'
#' It reads input of tree-level information (SPECIE, DBH, HT), stand-level (dominant specie,
#' zone), and one taper component (di or hi) to calculate the other element (hi or di) based on
#' a selection of available taper equation models.
#
#' @param SPECIES Dominant specie (1:Rauli, 2:Roble, 3:Coigue, 4:Others or Mixed)
#' @param ZONE Growth zone of the corresponding stand
#' @param DBH Diameter at breast height (cm) of tree
#' @param HT Total tree height (m)
#' @param di Stem diameter (cm) at given stem height hi (m)
#' @param hi Stem height (m) at given stem diameter di (cm)
#' @param model Type of selected taper model (1:zone specific, 2:all zones) (default = 2)
#'
#' @return The missing component (di or hi) from the requested taper model. Or it returns the complete
#' taper on vectors d and h (evaluated each centimeter up to total height).
#'
#' @author
#' S.A. Gezan and P. Moreno
#'
#' @references
#' Gezan, S.A. and Moreno, P. (2000). Informe Procedimientos y Resultados Modelos de Ahusamiento.
#' Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example 1: Unknown diameter inside bark stem diameter
#' (di <- get_taper(SPECIES=1, ZONE=2, DBH=8.1, HT=10.2, hi=4.3, model=1)$di)
#'
#' # Example 2: Unknown stem height
#' (hi <- get_taper(SPECIES=1, ZONE=2, DBH=8.1, HT=10.2, di=9.8, model=1)$hi)
#'
#' # Example 3: Plotting complete taper model for a given tree
#' tree <- get_taper(SPECIES=1, ZONE=2, DBH=8.1, HT=10.2, di=9.8, model=1)
#' plot(tree$h,tree$d, ylab='di (cm)', xlab='hi (m)', type='l', col='blue')

get_taper <- function(SPECIES = NA, ZONE=NA, DBH=NA, HT=NA, di=NA, hi=NA, model=2){

  if (is.na(model) || model==2){
    ZONE='Todas'
  }

  #Get parameters for that tree
  #taper_params_tree <- taper_params %>% filter(SPECIES. == SPECIES & ZONE. == ZONE)
  taper_params_tree <- taper_params[(taper_params$SPECIES.==SPECIES & taper_params$ZONE.==ZONE),]

  # This SPECIES. is with a dot. While the rest is without . is this confusing or can it cause errors?
  # Assigning to variables. It can be done with taper_params_tree$b0 always, but it is too long
  b0 <- taper_params_tree$b0
  b1 <- taper_params_tree$b1
  b2 <- taper_params_tree$b2
  b3 <- taper_params_tree$b3
  b4 <- taper_params_tree$b4
  b5 <- taper_params_tree$b5
  b6 <- taper_params_tree$b6
  b7 <- taper_params_tree$b7

  # Get the complete tree profile (in increments of 1 cm)
  h <- seq(0.01,HT,by=0.01)

  if(taper_params_tree$M == 4){
    x <- (HT-h)/(HT-1.3)
    y <- b0*x^(1.5) + b1*(x^(1.5) - x^3)*DBH + b2*(x^(1.5) - x^3)*HT + b3*(x^(1.5) - x^32)*DBH*HT + b4*(x^(1.5) - x^32)*HT^0.5 + b5*(x^(1.5) - x^40)*HT^2
    d <- DBH * sqrt(y)
  } else if (taper_params_tree$M == 5){
    z <- h/HT
    x <- (1-z^0.5) / (1-0.2^0.5)
    y <- b0 + b1*log(DBH) + b2*DBH + b3*log(x)*(z^2) + b4*log(x)*log(z) + b5*log(x)*(z^0.5) + b6*log(x)*exp(z) + b7*log(x)*(DBH/HT)
    d <- exp(y)
  } else if (taper_params_tree$M == 6){
    z <- h/HT
    x <- (1-z^0.5) / (1-0.2^0.5)
    y <- b0 + b1*log(DBH) + b2*log(x)*(z^2) + b3*log(x)*log(z) + b4*log(x)*(DBH/HT)
    d <- exp(y)
  }

  #Assigning last diameter to 0
  d[length(d)] <- 0

  # Getting missing parameter
  if (is.na(di)) {              # If di is missing
    posit <- which(abs(h - hi) == min(abs(h - hi)))
    di <- d[posit]
  } else if (is.na(hi)) {       # If hi is missing
    posit <- which(abs(d - di) == min(abs(d - di)))
    hi <- h[posit]
  }
  return(list(di=di,hi=hi,d=d,h=h))
}

# Note: - Need to check that di and hi are not illogical from the DBH and HT of the tree
# Needs to be corrected, taper_params has M 6, 5, no model 4, but required to be checked and specified by user!
