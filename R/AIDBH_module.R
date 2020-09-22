#' Projects individual annual increment in DBH for Nothofagus from different parameters
#'
#' Evaluates input variables to project individual annual increment in DBH for Nothofagus
#' based on tree-level and stand-level parameters. Projections are based in 1 year increments.
#'
#' @param SP Specie (1:Rauli, 2:Roble, 3:Coihue, 4:Other)
#' @param DBH Diameter at breast height (cm) of the tree
#' @param A Breast height age (year) of the tree
#' @param SS Sociological status (1:dominant, 2:codominant, 3:intermediate, 4:suppressed)
#' @param SSCAL Calculated sociological status
#' @param BALc Basal area of larger trees for Nothofagus (m2/ha)
#' @param ZONE Growth zone of the plot (1, 2, 3, 4)
#' @param AD Dominant age (years) of the plot
#' @param SDI Stand Density Index (trees/ha) of the plot
#' @param model Model to use for annual DBH increment (1, 2, 3, 4) (default = 1)
#'
#' @return Individual annual increment in DBH (AIDBH, mm) (for projection)
#'
#' @author
#' P. Moreno and S.A. Gezan
#'
#' @references
#' Moreno et al. (2017). Individual-tree diameter growth models for mixed Nothofagus
#' second growth forests in southern Chile. Forests 8(12), 506.
#'
#' @examples
#' # Example 1: Individual annual increment in DBH (AIDBH, mm) with model 1
#' Gest <- AIDBH_module(BALc=2, SDI=800, DBH=15, A=30, SS=1, model=1)
#' Gest
#'
#' # Example 2: Individual annual increment in DBH (AIDBH, mm) with model 2
#' Gest <- AIDBH_module(BALc=2, SDI=800, DBH=15, A=30, SS=1, AD=50, SSCAL=0.5, model=2)
#' Gest
#'
#' # Example 3: Individual annual increment in DBH (AIDBH, mm) with model 3
#' Gest <- AIDBH_module(SP=3, ZONE=2, BALc=2, SDI=800, DBH=15, A=30, SS=1, model=3)
#' Gest
#'
#' # Example 4: Individual annual increment in DBH (AIDBH, mm) with model 4
#' Gest <- AIDBH_module(SP=3, ZONE=2, BALc=2, SDI=800, DBH=15, A=30, SS=1, AD=50, SSCAL=0.5, model=4)
#' Gest

AIDBH_module <- function(SP=NA, DBH=NA, A=NA, SS=NA, SSCAL=NA, BALc=NA, ZONE=NA, AD=NA, SDI=NA, model=1){

  # Model 1 (CVselection): log(AIDBH) = b0 + b1*(BALc+10) + b2*SDI + b3*log(DBH) + b4*log(A) + b5*Ss
  betas1<-c(2.4097492437, -0.0070621315,  0.0002745541,  0.9045958975, -1.1380841801, -0.1335811948) # b0,b1,b2,b3,b4,b5

  # Model 2 (LASSO): log(AIDBH) = b0 + b1*Ss + b2*(BALc+10) + b3*SDI + b4*log(AD) + b5*log(DBH) + b6*log(A) + b7*(1/sqrt(PSCAL+10))
  betas2<-c(-0.3097512752, -0.1325217069, -0.0072549433, 0.0002215201, -0.0431837105, 0.7892183178, -1.0732148830, 9.7915258993) # b0,b1,b2,b3,b4,b5

  # Model 3 (CVselection with SPZONE): log(AIDBH) = b0 + b1*SPZONE12+ b2*SPZONE14+ b3*SPZONE21+ b4*SPZONE22+ b5*SPZONE23+ b6*SPZONE24+ b7*SPZONE31
  #           + b8*SPZONE32+ b9*SPZONE33+ b10*SPZONE34 + b11*(BALc+10) + b12*SDI + b13*log(DBH) + b14*log(A) + b15*Ss
  betas3<-c(2.65948954,0.20206391,0.35908797,-0.14256299,-0.09766417,0.15339191,-0.01030241,0.22600513,0.23241300,0.23972108,0.16592525,-0.00697406,0.00008419,0.91663284,-1.17524242,-0.14314288)

  # Model 4 (LASSO with SPZONE): log(AIDBH) = b0 + b1*Ss + b2*(BALc+10) + b3*SDI + b4*log(AD) + b5*log(DBH) + b6*log(A) + b7*(1/sqrt(PSCAL+10)
  # + b8*SPZONA14+ b9*SPZONA21+ b10*SPZONA22+ b11*SPZONA24+ b12*SPZONA32+ b13*SPZONA33)
  betas4<-c(-1.38682935,-0.12933416,-0.00711666,0.00008700,-0.02249610,0.76994689,-1.09622094,14.18970247,0.06790100,-0.25069820,-0.17015824,-0.19551723,0.03455293,0.05921630) # b0,b1,b2,b3,b4,b5

  ### Projection

  if (model==1){       ## Requires: BALc, SDI, DBH, A, SS
    if (is.na(BALc) || is.na(SDI) || is.na(DBH) || is.na(A) || is.na(SS)) {
      warning('Missing input information for AIDBH module evaluation', call.=FALSE)
    } else {
    est <- betas1[1] + betas1[2]*(BALc+10) + betas1[3]*SDI + betas1[4]*log(DBH) + betas1[5]*log(A) + betas1[6]*SS
    }
  }

  else if (model==2){  ## Requires: SS, BALc, SDI, AD, DBH, SSCAL, A
    if (is.na(SS) || is.na(BALc) || is.na(SDI) || is.na(AD) || is.na(DBH) || is.na(SSCAL) || is.na(A)) {
      warning('Missing input information for AIDBH module evaluation', call.=FALSE)
    } else {
    est <- betas2[1] + betas2[2]*SS + betas2[3]*(BALc+10) + betas2[4]*SDI + betas2[5]*log(AD) + betas2[6]*log(DBH) +
           betas2[7]*log(A) + betas2[8]*(1/sqrt(SSCAL+10))
    }
  }

  else if (model==3){  ## Requires: BALc, SDI, DBH, A, SS, and SPZONE 12, 14, 21, 22, 23, 24, 31, 32, 33, 34
    if (is.na(BALc) || is.na(SDI) || is.na(A) || is.na(SS) || is.na(DBH) || is.na(SP) || is.na(ZONE)) {
      warning('Missing input information for AIDBH module evaluation', call.=FALSE)
    } else {

    SPZONE12<-0; SPZONE14<-0; SPZONE21<-0; SPZONE22<-0; SPZONE23<-0
    SPZONE24<-0; SPZONE31<-0; SPZONE32<-0; SPZONE33<-0; SPZONE34<-0
    cod <- SP*10+ZONE
    if (cod==12) {SPZONE12=1}
    else if (cod==14){SPZONE14<-1}
    else if (cod==21){SPZONE21<-1}
    else if (cod==22){SPZONE22<-1}
    else if (cod==23){SPZONE23<-1}
    else if (cod==24){SPZONE24<-1}
    else if (cod==31){SPZONE31<-1}
    else if (cod==32){SPZONE32<-1}
    else if (cod==33){SPZONE33<-1}
    else {SPZONE34<-1}
    est <- betas3[1] + betas3[2]*SPZONE12 + betas3[3]*SPZONE14 + betas3[4]*SPZONE21 + betas3[5]*SPZONE22 +
                       betas3[6]*SPZONE23 + betas3[7]*SPZONE24 + betas3[8]*SPZONE31 + betas3[9]*SPZONE32 +
                       betas3[10]*SPZONE33 + betas3[11]*SPZONE34 +
                       betas3[12]*(BALc+10) + betas3[13]*SDI + betas3[14]*log(DBH) + betas3[15]*log(A) + betas3[16]*SS
    }
  }

  else if (model==4){  ## Requires: SS, BALc, SDI, AD, DBH, A, SSCAL, and SPZONE 14, 21, 22, 24, 32, 33
    if (is.na(BALc) || is.na(SDI) || is.na(A) || is.na(SS) || is.na(DBH) || is.na(AD) || is.na(SSCAL) || is.na(SP) || is.na(ZONE)) {
      warning('Missing input information for AIDBH module evaluation', call.=FALSE)
    } else {

    SPZONE14<-0; SPZONE21<-0; SPZONE22<-0; SPZONE24<-0; SPZONE32<-0; SPZONE33<-0
    cod<-SP*10+ZONE
    if (cod==14) {SPZONE14<-1}
    else if (cod==21){SPZONE21<-1}
    else if (cod==22){SPZONE22<-1}
    else if (cod==24){SPZONE24<-1}
    else if (cod==32){SPZONE32<-1}
    else if (cod==33){SPZONE33<-1}
    est <- betas4[1] +  betas4[2]*SS + betas4[3]*(BALc+10) + betas4[4]*SDI + betas4[5]*log(AD) +
                        betas4[6]*log(DBH) + betas4[7]*log(A) + betas4[8]*(1/(sqrt(SSCAL+10))) +
                        betas4[9]*SPZONE14 + betas4[10]*SPZONE21 + betas4[11]*SPZONE22 + betas4[12]*SPZONE24 +
                        betas4[13]*SPZONE32 + betas4[14]*SPZONE33
    }
  }

  AIDBH <- exp(est)
  return(AIDBH)

}
