#' Simulates plot at the whole stand-level for the next period
#'
#' It simulates plot level growth, mortality (and recruitment) of a given stand
#' requiring stand-level parameters coming from core_module.
#' Simulations are done using stand-level models starting from intial age (AD0)
#' until final age (ADF) in increments of 1 year.
#' Note that recruitment is currently not considered.
#'
#' @param core.stand An object originated from core_module
#'
#' @return A series of elements and parameters with updated tables. Outputs are the same as \emph{input_module}
#'         but these are updated for simulations. The main table is \emph{sp.table} that contains a stand table
#'         for the new simulations for the given parameters. Also, all other information is updated for
#'         time AF.
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example: Input stand-level data to predict volume using stand-level functions to next year
#' BA <- c(36.5, 2.8, 1.6, 2.4)
#' N <- c(464, 23, 16, 48)
#' input <- input_module(ZONE=1, AD=28, HD=18.5, AF=29, N=N, BA=BA,
#'                       type='stand', ddiam=FALSE, V_model=1)
#' input$sp.table  # Current stand at AD (time 0)
#' current.plot <- stand_simulator(core.stand=input)
#' current.plot$sp.table # Stand at AF=AD+1

stand_simulator <- function(core.stand=NULL){

  #Error with dominant species
  if (core.stand$DOM.SP == 4){
    stop('The stand is not dominated by Nothofagus. The simulator only works with Nothofagus dominated forests.')
  }

  core.stand$sp.table$PSP.NHA<-core.stand$sp.table[,2]/sum(core.stand$sp.table[1:3,2])
  core.stand$sp.table$PSP.BA<-core.stand$sp.table[,3]/sum(core.stand$sp.table[1:3,3])
  #print(core.stand$sp.table)

  SI <- core.stand$SI
  DOM.SP <- core.stand$DOM.SP
  ZONE <- core.stand$ZONE
  HD0 <- core.stand$HD

  #NHAN0 <- core.stand$sp.table$N[1:3] %>% sum(na.rm = TRUE)
  NHAN0 <- sum(core.stand$sp.table$N[1:3], na.rm=TRUE)
  NHA0 <- core.stand$sp.table$N[5]
  NHA990 <- NHA0-NHAN0

  QD0 <- core.stand$sp.table$QD[5]
  BAN0 <- sum(core.stand$sp.table$BA[1:3], na.rm=TRUE)
  #BAN0 <- core.stand$sp.table$BA[1:3] %>% sum(na.rm = TRUE)
  BA990 <- core.stand$sp.table$BA[4]
  BA0 <- BAN0+BA990

  PBAN0 <- core.stand$PBAN
  PNHAN0 <- core.stand$PNHAN

  PBAN1 <- PBAN0   # Next year is same as current year - both updated after simulation
  PNHAN1 <- PNHAN0 # Next year is same as current year

    # Simulates from y=AD to y=AD+1
    y<-core.stand$AD
    NHA1 <- NHAmodule(NHA0=NHA0, QD0=QD0, DOM.SP=DOM.SP, model=core.stand$NHA_model)   #Estimates new number of trees
    PBAN1 <- PBAN0   # Next year is same as current year (this needs future work)
    BAN1 <- BANmodule(BAN0=BAN0, AD0=y, SI=SI, NHA0=NHA0, NHA1=NHA1, PBAN0=PBAN0, PBAN1=PBAN1, projection=TRUE)$BAN1   #projects new basal area (needs to change)

    if (PNHAN0 == 1) {
      PNHAN1 <- PNHAN0
    } else {
      # Update PNHAN1 - based on the differences on time, otherwise, it is a very drastic change
      Temp_PNHAN0 <- exp(-7.13684+10.29084*PBAN0-0.01404*(y))/(1+exp(-7.13684+10.29084*PBAN0-0.01404*(y)))
      Temp_PNHAN1 <- exp(-7.13684+10.29084*PBAN1-0.01404*(y+1))/(1+exp(-7.13684+10.29084*PBAN1-0.01404*(y+1)))
      PNHAN1 <- PNHAN0 + (Temp_PNHAN0-Temp_PNHAN1) # Important change
    }
    if (PBAN0 == 1) {
      BA991 <- 0
    } else {
      BA991 <- BA99module(BA990=BA990, AD0=y, PNHAN0=PNHAN0, PNHAN1=PNHAN1, PBAN0=PBAN0, PBAN1=PBAN1, projection=TRUE)$BA991   #projects new basal area (needs to change)
    }
    BA1 <- BAN1 + BA991 # Calculates total new Basal Area
    PBAN1 <- BAN1/BA1   # Updates PBAN1
    NHAN1 <- NHA1*PNHAN1
    NHA991 <- NHA1*(1-PNHAN1)

    QD1 <- get_stand(BA=BA1, N=NHA1)   #New quadratic diameter
    HD1 <- get_site(DOM.SP=DOM.SP, ZONE=ZONE, SI=SI, AD=y)   #New dominant height

  #sp.table is the final stand values (for each year)
  sp.table<-matrix(data=0,nrow=5,ncol=4)
  sp.table[,1]<-c(seq(1:4),0)
  sp.table[1,2]<-core.stand$sp.table$PSP.NHA[1]*NHAN1
  sp.table[2,2]<-core.stand$sp.table$PSP.NHA[2]*NHAN1
  sp.table[3,2]<-core.stand$sp.table$PSP.NHA[3]*NHAN1
  sp.table[4,2]<-NHA991
  sp.table[5,2]<-NHA1

  sp.table[1,3]<-core.stand$sp.table$PSP.BA[1]*BAN1
  sp.table[2,3]<-core.stand$sp.table$PSP.BA[2]*BAN1
  sp.table[3,3]<-core.stand$sp.table$PSP.BA[3]*BAN1
  sp.table[4,3]<-BA991
  sp.table[5,3]<-BA1

  sp.table[1,4]<-get_stand(BA=sp.table[1,3], N=sp.table[1,2])
  sp.table[2,4]<-get_stand(BA=sp.table[2,3], N=sp.table[2,2])
  sp.table[3,4]<-get_stand(BA=sp.table[3,3], N=sp.table[3,2])
  sp.table[4,4]<-get_stand(BA=sp.table[4,3], N=sp.table[4,2])
  sp.table[5,4]<-get_stand(BA=sp.table[5,3], N=sp.table[5,2])

  sp.table<-data.frame(sp.table)
  colnames(sp.table)<-c('SPECIE','N','BA','QD')

  output <- input_module(ZONE=core.stand$ZONE,
                         AD=y+1, AF=core.stand$AF, HD=HD1, SI=core.stand$SI,
                         N=sp.table$N[1:4], BA=sp.table$BA[1:4],
                         area=core.stand$area, type=core.stand$type, ddiam=core.stand$ddiam, comptype=core.stand$comptype,
                         NHA_model=core.stand$NHA_model, V_model=core.stand$V_model, T_model=core.stand$T_model,
                         IADBH_model=core.stand$IADBH_model,
                         ATHIN=core.stand$ATHIN, BARp=core.stand$BARp, QD_ba=core.stand$QD_ba,
                         tree.list=NA)

  return(output)

}
