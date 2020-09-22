#' Estimates individual-tree volume inside bark based on taper equations
#'
#' It reads input of tree information (specie, DBH, HT), stand information (SPECIES, ZONE),
#' and of specification of restrictions (dmin, blength), to estimate inside bark volume (m3)
#' for the given tree based on different fitted taper equation models. Tree volume is calculated
#' starting on stump (0.3 m) until one of the following is reached: 1) stem diameter limit (dmin),
#' or 2) bole lenght (blength)
#'
#' @param ZONE Growth zone (1, 2, 3, 4) of the stand
#' @param SPECIES Species of the given tree (1:Rauli, 2:Roble, 3:Coigue)
#' @param DBH Diameter at breast height (cm)
#' @param HT Total tree height (m)
#' @param dmin Minimum stem diameter inside bark (cm)
#' @param blength Bole length (m)
#' @param stump Length of stump to discount (default = 0.3 m)
#' @param model Model selected for taper equation (1:zone specific, 2:all zones) (default = 2)
#'
#' @return Volume of stem section without bark (m3) based on the restrictions provided
#'
#' @author
#' S.A. Gezan & P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' Gezan, S.A. and Moreno, P. (2000b). INFORME MODELOS DE VOLUMEN.
#' Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example 1: Calculates tree volume for diameter limit of 5 cm (with stump of 0.3 m)
#' Vmodule_individual(SPECIES=1, ZONE=1, DBH=22.1, HT=18.2, dmin=5, model=1)
#'
#' # Example 2: Calculates tree volume for a bole length of 6 m (with a stump of 0.3 m)
#' Vmodule_individual(SPECIES=1, ZONE=1, DBH=22.1, HT=18.2, blength=6, model=1)
#'
#' # Example 3: Two ways to calculate total tree volume (with a stump of 0.3 m)
#' Vmodule_individual(SPECIES=1, ZONE=1, DBH=22.1, HT=18.2, dmin=0, model=1)
#' Vmodule_individual(SPECIES=1, ZONE=1, DBH=22.1, HT=18.2, blength=18.2, model=1)
#'
#' # Example 4: Calculates total tree volume without discounting for stump
#' Vmodule_individual(SPECIES=1, ZONE=1, DBH=22.1, HT=18.2, dmin=0, stump=0, model=1)

Vmodule_individual <- function(SPECIES=NA, ZONE=NA, DBH=NA, HT=NA, dmin=NA,
                               blength=NA, stump=0.3, model=2){

  if (is.na(blength) & is.na(dmin)){
    stop('Minimum diameter or bole length need to be provided.')
  }
  if (!is.na(blength) & !is.na(dmin)){
    stop('Minimum diameter or bole length are both provided.')
  }
  # dmin provided
  if (is.na(blength) & !is.na(dmin)){
    blength<-get_taper(SPECIES=SPECIES, ZONE=ZONE, DBH=DBH, HT=HT, di=dmin, model=model)$hi
  }

  incr<-0.01  # default increment in h from get_taper
  tree.profile<-get_taper(SPECIES=SPECIES, ZONE=ZONE, DBH=DBH, HT=HT, hi=blength, model=model)
  d<-tree.profile$d
  h<-tree.profile$h
  ba0<-pi*((d/100)^2)/4
  ba1<-c(ba0[2:length(ba0)],0)
  vsection<-incr*(ba0+ba1)/2

  s.start<-which.min(abs(h-stump))
  s.fin<-which.min(abs(h-blength))
  vtree<-sum(vsection[s.start:s.fin])

  return(vtree)
}

# Note
# - If we need a different product stump is the lower height and blength is the upper height
# - Remember it is di without bark, hence at hi=1.3, di<>DBH.
