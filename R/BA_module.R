#' Predicts and/or projects stand Basal Area (BA) from stand-level parameters
#'
#' Evaluates stand-level input variables to predict and/or project basal area.
#' Projections are based in 1 year increments.
#'
#' @param AD0 Dominant age (years) of the stand at current time
#' @param HD0 Dominant height (m) of the stand at current time
#' @param N0 Number of trees (trees/ha) of the stand at current time
#' @param BA0 Basal area (m2/ha) of the stand at current time (required for projections)
#' @param HD1 Dominant height (m) of the stand at future projection time
#' (required for projections)
#' @param N1 Number of trees (trees/ha) of the stand at future projection time
#' (required for projections)
#' @param model Model to use for estimation (1:non-linear fit, 2:linear fit) (default = 1)
#' @param projection if TRUE projection from BA0 is executed for a 1 year increment
#'
#' @return Basal area (BA0, m2/ha) for the current age (for prediction)
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
#' # Example 1: Predicts basal area at a certain age AD
#' BAest <- BAmodule(AD0=20, HD0=17.2, N0=2730, model=1, projection=FALSE)
#' BAest$BA0
#'
#' # Example 2: Projects basal area for different models
#' (BAest <- BAmodule(AD0=20, HD0=17.2, N0=2730, BA0=33.11, HD1=19.1,
#'                    N1=2610, model=1, projection=TRUE)$BA1)
#' (BAest <- BAmodule(AD0=20, HD0=17.2, N0=2730, BA0=33.11, HD1=19.1,
#'                    N1=2610, model=2, projection=TRUE)$BA1)

BAmodule <- function(AD0=NA, HD0=NA, N0=NA, BA0=NA, HD1=NA, N1=NA, model=1, projection=FALSE){

  # Model 1 (non-linear): BA = exp(b0)*AD^b1*HD^b2*NHA^b3
  bm1<-c(-3.36953,0.45753,0.75578,0.41176) # b0,b1,b2,b3
  # Model 2 (linear): log(AB) = b0 + b1*log.AD + b2*log.HD + b3*log.NHA
  bm2<-c(-3.90894,0.43824,0.87593,0.43949) # b0,b1,b2,b3
  if (model==1){
    bm<-bm1
  }
  if (model==2){
    bm<-bm2
  }

  if (projection==FALSE){
    # Prediction
    BA0<-exp(bm[1])*(AD0^bm[2])*(HD0^bm[3])*(N0^bm[4])
    BA1<-NA
  }
  if (projection==TRUE){
    # Projection
    #SI <- get_site(dom_sp=1, zone=1, HD=HD0, AD=AD0)
    #HDA <- get_site(dom_sp=1, zone=1, SI=SI, AD=AD0-0.5)
    #HDB <- get_site(dom_sp=1, zone=1, SI=SI, AD=AD0+0.5)
    #DerHD<-(HDB-HDA) # Derivative HD
    #BA1<-BA0*(1+bm[2]/AD0)  # Needs to be adjusted for other derivatives
    #BA1<-BA0*(1+(bm[2]/AD0)+(bm[3]/HD0)*DerHD)+(bm[4]/N0)*DerN)  # Needs to be adjusted for other derivatives
    #derHD es DerN. No sabemos la derivada con el numero de arbole spor hectarea
    #la derivad apodria ser 98 y regreso al modelo y uso un valor de 98 y da una nueva area basal.
    #medio iterativa

    #BA1<-BA0exp(bm[2]*log((AD0+1)/AD0))  # Needs to be adjusted for other derivatives
    BA1<-BA0*exp(bm[2]*log((AD0+1)/AD0)+bm[3]*log(HD1/HD0)+bm[4]*log(N1/N0))
  }

  return(list(BA0=BA0,BA1=BA1))
}

# Note: BA1 model needs to incorporate other derivatives in the future
