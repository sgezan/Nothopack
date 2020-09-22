#' Generates the diameter distribution for a given Nothofagus stand to each species/cohort
#'
#' It generates the diameter distribution for a given stand for each of the
#' species/cohorts using the method of parameter recovery. The diameter classes are based on
#' 5 cm intervals (by default). Data can be read from \emph{input_module}.
#'
#' @param sp.table table with stand-level information by specie and total. Columns are
#'            SPECIES, N, BA, QD, and table is grouped by SPECIES (1:Rauli, 2:Roble, 3:Coigue, 4:Other, 0:Total).
#' @param HD Dominant height (m) of dominant specie in the current stand
#' @param DOM.SP Dominant specie (1: Rauli, 2: Roble, 3: Coigue, 4:Mixed)
#' @param ZONE Growth zone (1, 2, 3, 4)
#' @param class Diameter class width (cm) (default = 5)
#'

#' @return A 3-dimension matrix of values for each specie/cohort by diameter classes (on specified
#' increments) by SPECIES and then by all relevant columns.
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example: Generation of distribution from stand-level input data
#' BA <- c(36.5, 2.8, 0, 2.4)
#' N <- c(464, 23, 0, 48)
#' plot <- input_module(type='stand', ZONE=2, AD=28, AF=28, HD=23.5, N=N, BA=BA)
#' Dd <- diam_dist(sp.table=plot$sp.table, HD=plot$HD,
#'               DOM.SP=plot$DOM.SP, ZONE=plot$ZONE)
#' Dd[5,,]  # Total diameter distribution
#' # Ploting distribution for all species
#' barplot(as.matrix(Dd[5,,5]), main='Diameter Distribution all species',
#'         xlab='DBH Class', beside=TRUE, col=4)

diam_dist <- function(sp.table=NA, HD=NA, DOM.SP=NA, ZONE=NA, class=5){

  #print(sp.table)
  vNHA <- sp.table$N
  vBA <- sp.table$BA
  vQD <- sp.table$QD

  BA <- vBA[5]
  NHA <- vNHA[5]
  QD <- vQD[5]

  NHAN <- sum(vNHA[1:3])   # Total number of trees Nothofagus
  PNHAN <- NHAN/NHA        # Proportion trees Nothofagus

  A <- 5  # minimum diamater for any distribution

  # Rauli parameters, sp=1
  b0_B_1 <- -4.85344986
  b1_B_1 <-  1.03816837
  b2_B_1 <- -0.01755607
  b3_B_1 <- 11.88124162
  b0_C_1 <-  3.41813315
  b1_C_1 <-  0.46054858
  b2_C_1 <- -0.42701842

  # Roble parameters, sp=2
  b0_B_2 <- -8.83509885
  b1_B_2 <-  1.06360951
  b2_B_2 <-  0.99954695
  b0_C_2 <-  3.57525492
  b1_C_2 <-  0.39042726
  b2_C_2 <- -0.36299573
  b3_C_2 <- -4.83178811

  #Coigue parameters, sp=3
  b0_B_3 <- -6.57263098
  b1_B_3 <-  1.09071507
  b0_C_3 <-  5.72487482
  b1_C_3 <-  0.64094474
  b2_C_3 <- -0.64908497
  b3_C_3 <- -9.04444521

  # Other sps parameters, sp=4
  b0_B_4 <- -2.06868340
  b1_B_4 <-  0.68162016
  b2_B_4 <- -0.00029245
  b0_C_4 <-  2.15124820
  b1_C_4 <-  0.29765454
  b2_C_4 <- -0.22722466

  # Rauli
  B1 <- b0_B_1 + b1_B_1*vQD[1] + b2_B_1*PNHAN + b3_B_1/BA
  C1 <- b0_C_1 + b1_C_1*B1 + b2_C_1*vQD[1]
  # Roble
  RS <- 100 * sqrt(10000/NHA)/HD
  B2 <- b0_B_2 + b1_B_2*vQD[2] + b2_B_2*log(RS)
  C2 <- b0_C_2 + b1_C_2*B2 + b2_C_2*vQD[2] + b3_C_2/QD
  # Coigue
  B3 <- b0_B_3 + b1_B_3*vQD[3]
  C3 <- b0_C_3 + b1_C_3*B3 + b2_C_3*vQD[3] + b3_C_3/QD
  # Others
  B4 <- b0_B_4 + b1_B_4*vQD[4] + b2_B_4*NHAN
  C4 <- b0_C_4 + b1_C_4*B4 + b2_C_4*vQD[4]

  # Calculation of the probability of NHA per diameter class, from 5 to 80 cm

  DBH_LL <- matrix(data=0,nrow=0,ncol=1)
  DBH_UL <- matrix(data=0,nrow=0,ncol=1)
  Prob1 <- matrix(data=0,nrow=0,ncol=1)
  Prob2 <- matrix(data=0,nrow=0,ncol=1)
  Prob3 <- matrix(data=0,nrow=0,ncol=1)
  Prob4 <- matrix(data=0,nrow=0,ncol=1)
  Dclass <- matrix(data=0,nrow=0,ncol=1)
  BAclass <- matrix(data=0,nrow=0,ncol=1)
  Vclass <- matrix(data=0,nrow=0,ncol=1)
  Hclass1 <- matrix(data=0,nrow=0,ncol=1)
  Hclass2 <- matrix(data=0,nrow=0,ncol=1)
  Hclass3 <- matrix(data=0,nrow=0,ncol=1)
  Hclass4 <- matrix(data=0,nrow=0,ncol=1)
  Hclass5 <- matrix(data=0,nrow=0,ncol=1)

  # Recover of diameter distribution parameters by species

  maxD <- 95/class
  diam <- seq(from=5,to=95,by=class)   # Diameter classes

  for (j in 1:(length(diam)-1)){
    DBH_LL[j] <- diam[j]                 # cm
    DBH_UL[j] <- diam[j+1]               # cm
    Dclass[j] <- (diam[j]+diam[j+1])/2   # cm
    BAclass[j] <- (pi/4)*((Dclass[j])^2) # cm2

    #Hclass1[j] <- height_param(HD=HD, QD=QD, DBH=Dclass[j], dom_sp=1, zone=zone)
    #Hclass2[j] <- height_param(HD=HD, QD=QD, DBH=Dclass[j], dom_sp=2, zone=zone)
    #Hclass3[j] <- height_param(HD=HD, QD=QD, DBH=Dclass[j], dom_sp=3, zone=zone)
    #Hclass4[j] <- height_param(HD=HD, QD=QD, DBH=Dclass[j], dom_sp=DOM.SP, zone=zone)

    Hclass1[j] <- height_param(DOM.SP=1, ZONE=ZONE, HD=HD, QD=QD, DBH=Dclass[j])
    Hclass2[j] <- height_param(DOM.SP=2, ZONE=ZONE, HD=HD, QD=QD, DBH=Dclass[j])
    Hclass3[j] <- height_param(DOM.SP=3, ZONE=ZONE, HD=HD, QD=QD, DBH=Dclass[j])
    Hclass4[j] <- height_param(DOM.SP=DOM.SP, ZONE=ZONE, HD=HD, QD=QD, DBH=Dclass[j])

    #height_param(DOM.SP=2, ZONE=2, HD=15, QD=12, DBH=24)

    if (Hclass1[j]<1.3) { Hclass1[j]<-1.3 }
    if (Hclass2[j]<1.3) { Hclass2[j]<-1.3 }
    if (Hclass3[j]<1.3) { Hclass3[j]<-1.3 }
    if (Hclass4[j]<1.3) { Hclass4[j]<-1.3 }

    if (vNHA[1]==0) {
      Prob1[j]=0
      #Hclass1[j]=0
    } else { Prob1[j] <- exp(-((diam[j] - A)/B1)^C1) - exp(-((diam[j+1] - A)/B1)^C1) }

    if (vNHA[2]==0) {
      Prob2[j]=0
      #Hclass2[j]=0
    } else { Prob2[j] <- exp(-((diam[j] - A)/B2)^C2) - exp(-((diam[j+1] - A)/B2)^C2) }

    if (vNHA[3]==0) {
      Prob3[j]=0
      #Hclass3[j]=0
    } else { Prob3[j] <- exp(-((diam[j] - A)/B3)^C3) - exp(-((diam[j+1] - A)/B3)^C3) }

    if (vNHA[4]==0) {
      Prob4[j]=0
      #Hclass4[j]=0
    } else { Prob4[j] <- exp(-((diam[j] - A)/B4)^C4) - exp(-((diam[j+1] - A)/B4)^C4) }

    if (is.na(Prob1[j])) {Prob1[j] = 0}
    if (is.na(Prob2[j])) {Prob2[j] = 0}
    if (is.na(Prob3[j])) {Prob3[j] = 0}
    if (is.na(Prob4[j])) {Prob4[j] = 0}

  }

  # Adjust probabilities for truncation with p = 0.05 (p% on the upper tail is truncated)
  p <- 0.05
  CProb1 <- cumsum(Prob1)/(1-p)
  CProb2 <- cumsum(Prob2)/(1-p)
  CProb3 <- cumsum(Prob3)/(1-p)
  CProb4 <- cumsum(Prob4)/(1-p)

  Prob1[1] <- CProb1[1]
  Prob2[1] <- CProb2[1]
  Prob3[1] <- CProb3[1]
  Prob4[1] <- CProb4[1]

  for (j in 2:(length(Prob1))){
    if (CProb1[j]>1) {CProb1[j]=1}
    if (CProb2[j]>1) {CProb2[j]=1}
    if (CProb3[j]>1) {CProb3[j]=1}
    if (CProb4[j]>1) {CProb4[j]=1}

    Prob1[j] <- CProb1[j] - CProb1[j-1]
    Prob2[j] <- CProb2[j] - CProb2[j-1]
    Prob3[j] <- CProb3[j] - CProb3[j-1]
    Prob4[j] <- CProb4[j] - CProb4[j-1]

  }

  N.sp1 <- round(vNHA[1]*Prob1,3) # only 1 decimal?
  N.sp2 <- round(vNHA[2]*Prob2,3)
  N.sp3 <- round(vNHA[3]*Prob3,3)
  N.sp4 <- round(vNHA[4]*Prob4,3)
  N.total <- N.sp1 + N.sp2 + N.sp3 + N.sp4

  # Adjusting N.total with NHA
  if (vNHA[1]==0) { K1=0
  } else {
    K1 <- sum(N.sp1)/vNHA[1]
    N.sp1 <- round(N.sp1/K1,3)
  }
  if (vNHA[2]==0) { K2=0
  } else {
    K2 <- sum(N.sp2)/vNHA[2]
    N.sp2 <- round(N.sp2/K2,3)
  }
  if (vNHA[3]==0) { K3=0
  } else {
    K3 <- sum(N.sp3)/vNHA[3]
    N.sp3 <- round(N.sp3/K3,3)
  }
  if (vNHA[4]==0) { K4=0
  } else {
    K4 <- sum(N.sp4)/vNHA[4]
    N.sp4 <- round(N.sp4/K4,3)
  }

  BA.sp1 <- round(BAclass*N.sp1/(100^2),5)  # cm2/ha
  BA.sp2 <- round(BAclass*N.sp2/(100^2),5)  # cm2/ha
  BA.sp3 <- round(BAclass*N.sp3/(100^2),5)  # cm2/ha
  BA.sp4 <- round(BAclass*N.sp4/(100^2),5)  # cm2/ha
  BA.total <- BA.sp1 + BA.sp2 + BA.sp3 + BA.sp4

  # Adjusting BA.total with BA
  if (vBA[1]==0) { K1=0
  } else {
    K1 <- sum(BA.sp1)/vBA[1]
    BA.sp1 <- round(BA.sp1/K1,5)
  }
  if (vBA[2]==0) { K2=0
  } else {
    K2 <- sum(BA.sp2)/vBA[2]
    BA.sp2 <- round(BA.sp2/K2,5)
  }
  if (vBA[3]==0) { K3=0
  } else {
    K3 <- sum(BA.sp3)/vBA[3]
    BA.sp3 <- round(BA.sp3/K3,5)
  }
  if (vBA[4]==0) { K4=0
  } else {
    K4 <- sum(BA.sp4)/vBA[4]
    BA.sp4 <- round(BA.sp4/K4,5)
  }

  BA.total <- BA.sp1 + BA.sp2 + BA.sp3 + BA.sp4

  # Adjusting heights

  # Observed HD from estimated Hclass (for any specie using Hclass5)
  Hclass5 <- (Hclass1*N.sp1 + Hclass2*N.sp2 + Hclass3*N.sp3 + Hclass4*N.sp4)/(N.total)
  Hclass5 <- replace(Hclass5, is.na(Hclass5), 0)

  Temp.Prod <- Hclass5*N.total
  Temp.Sum <- 0
  Ncount <- 0
  for (j in (length(diam)-1):1) {
    Ncount <- Ncount + N.total[j]
    if (Ncount > 100) {
      Htemp <- (Temp.Sum+(Hclass5[j]*(100-(Ncount-N.total[j]))))/100
      break
    } else {
      Temp.Sum <- Temp.Sum + Temp.Prod[j]
    }
  }

  K <- HD/Htemp
  Hclass1 <- round(K*Hclass1,2)
  Hclass2 <- round(K*Hclass2,2)
  Hclass3 <- round(K*Hclass3,2)
  Hclass4 <- round(K*Hclass4,2)

  Hclass5 <- (Hclass1*N.sp1 + Hclass2*N.sp2 + Hclass3*N.sp3 + Hclass4*N.sp4)/(N.total)
  Hclass5 <- replace(Hclass5, is.na(Hclass5), 0)
  r_names<-c('DBH_ll','DBH_ul','D_class','H_class','N','BA')
  DDist<-array(data=NA, dim=c(5,(length(diam)-1),6),
               dimnames=list(c(1:5),c(1:(length(diam)-1)),r_names))
  DDist[1,,]<-cbind(DBH_LL,DBH_UL,Dclass,Hclass1,N.sp1,BA.sp1)  #1: Rauli
  DDist[2,,]<-cbind(DBH_LL,DBH_UL,Dclass,Hclass2,N.sp2,BA.sp2)  #2: Roble
  DDist[3,,]<-cbind(DBH_LL,DBH_UL,Dclass,Hclass3,N.sp3,BA.sp3)  #3: Coigue
  DDist[4,,]<-cbind(DBH_LL,DBH_UL,Dclass,Hclass4,N.sp4,BA.sp4)  #4: Others
  DDist[5,,]<-cbind(DBH_LL,DBH_UL,Dclass,Hclass5,N.total,BA.total)  #0: Total

  # Add-up to make H_Class Zero on DDist[1,,4]
  for (j in (1:nrow(DDist[1,,]))) {
    if (DDist[1,j,5]==0) { DDist[1,j,4]=0 }
    if (DDist[2,j,5]==0) { DDist[2,j,4]=0 }
    if (DDist[3,j,5]==0) { DDist[3,j,4]=0 }
    if (DDist[4,j,5]==0) { DDist[4,j,4]=0 }
    if (DDist[5,j,5]==0) { DDist[5,j,4]=0 }
  }

  return(stand.table=DDist)

}

# Note: - Need to check that PNHAN goes from 0 to 1, and not 0 to 100
# It will be good to set up a minimum numer of trees per class to make tail shorter...
#' dom_sp and zone should be specified
