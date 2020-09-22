#' Inventory processing by calculating stand-level parameter based on individual tree data
#'
#' Calculates all relevant stand-level parameters from an inventory plot for
#' each of the species and for all together (1:Rauli, 2:Roble, 3:Coigue, 4:Other, 0:All).
#' The input table should contain the colums of SPECIES (1:Rauli, 2:Roble, 3:Coigue, 4:Other)
#' together with their diameter at breast height (DBH, cm), total tree height (HT, m), and
#' expansion factor (FT = 10,000/area, ha/m2). If area is provided, then expansion factor is not
#' considered.
#'
#' @param plotdata Data frame with plot data including the columns: SPECIES (1, 2, 3, 4), DBH (cm) and HT (m)
#' @param area Area of plot (m2)
#' @param extra List of additional column numbers to process from plot data for summary on stand table
#'
#' @return A data frame with columns stand-table:
#'                                 SPECIES (1:Rauli, 2:Roble, 3:Coigue, 4:Others, 0:All),
#'                                 N (number of trees/ha), BA (basal area, m2/ha), and QD (quadratic diameter, cm).
#'         Also stand parameters:   DOM.SP (dominant specie), HD (dominant height, m),
#'                                 PBAN (proportion of basal area of Nothofagus),
#'                                 PNHAN (proportion of number of trees of Nothofagus), and
#'                                 BAind (basal area for each individual tree, m2)
#'         If an extra list is provided, then additional summary columns are added to the
#'         stand-table grouped by specie.
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example 1: Simple plot calculations
#' plotdata <- plot_example
#' head(plotdata,10)
#' plot.temp <- input_module(ZONE=1, AD=28, type='tree', area=500,
#'                           tree.list=plotdata, T_model=2)$tree.list
#' head(plot.temp,10)
#' plot.inv <- stand_parameters(plotdata=plot.temp, area=500)
#' plot.inv$DOM.SP
#' plot.inv$HD
#' plot.inv$sd
#'
#' # Example 2: Plot calculations extended to additional columns
#' plotdata <- plot_example
#' head(plotdata,10)
#' plot.temp <- input_module(ZONE=1, AD=28, type='tree', area=500,
#'                           tree.list=plotdata, T_model=2)$tree.list
#' head(plot.temp,10)
#' plot.extra <- stand_parameters(plotdata=plot.temp, area=500, extra=8)
#' plot.extra$sd

stand_parameters <- function(plotdata=NA, area=NA, extra=NA){

  if (sum(is.na(plotdata$FT))>0) {
    FT<-matrix(1,nrow=nrow(plotdata),ncol=1)
    plotdata$FT<-FT*(10000/area)
  }

  # Calculating individual basal area
  plotdata$BAind <- plotdata$FT * as.numeric(pi * plotdata$DBH^2/40000 )  # Units: m2

  # Adding a factor to complete
  plotdata$SPECIES = factor(plotdata$SPECIES, levels=c(1:4,0))

  # table with N, BA and QD. Adds 0 if no species is found
  # sd <- plotdata %>% group_by(SPECIES) %>%
  #   summarise(
  #     N = sum(FT, na.rm = TRUE),
  #     BA = sum(BAind, na.rm = TRUE),
  #     QD = get_stand(BA, N)
  #   ) %>%
  #   complete(SPECIES, fill = list(N = 0, BA = 0, QD = 0)) %>%
  #   as.data.frame()
  sd<-matrix(0,nrow=5,ncol=4)
  sd[,1]<-c(1,2,3,4,0)
  colnames(sd)<-c('SPECIES','N','BA','QD')
  sd<-as.data.frame(sd)

  sd[1,2] <- sum(plotdata$FT[plotdata$SPECIES==1], na.rm=TRUE)
  sd[2,2] <- sum(plotdata$FT[plotdata$SPECIES==2], na.rm=TRUE)
  sd[3,2] <- sum(plotdata$FT[plotdata$SPECIES==3], na.rm=TRUE)
  sd[4,2] <- sum(plotdata$FT[plotdata$SPECIES==4], na.rm=TRUE)
  sd[1,3] <- sum(plotdata$BAind[plotdata$SPECIES==1], na.rm=TRUE)
  sd[2,3] <- sum(plotdata$BAind[plotdata$SPECIES==2], na.rm=TRUE)
  sd[3,3] <- sum(plotdata$BAind[plotdata$SPECIES==3], na.rm=TRUE)
  sd[4,3] <- sum(plotdata$BAind[plotdata$SPECIES==4], na.rm=TRUE)
  sd[1,4] <- get_stand(sd[1,3], sd[1,2])
  sd[2,4] <- get_stand(sd[2,3], sd[2,2])
  sd[3,4] <- get_stand(sd[3,3], sd[3,2])
  sd[4,4] <- get_stand(sd[4,3], sd[4,2])

  # Adding total values
  sd[5,2:4] <- c(sum(sd$N, na.rm = TRUE),
                 sum(sd$BA, na.rm = TRUE),
                 get_stand(BA = sum(sd$BA, na.rm = TRUE), N = sum(sd$N, na.rm = TRUE)))

  #proportion values
  PBAN <- sum(sd$BA[1:3])/sd$BA[5]
  PNHAN <- sum(sd$N[1:3])/sd$N[5]
  DOM.SP<-get_domsp(BA=sd$BA[1:4])

  # Dominant Height - 100 trees with largest DBH
  # (this is for any of the SPECIES, not only dominant sp)
  # It uses the CF or FT for the original plot, not the updated one form simulations.
  if (sum(is.na(plotdata$HT)>0)) {
    warning('HD can not be calculated as some HT are missing', call.=FALSE)
    HD <- NA
  } else {
    CF <- 10000/area  # Correction factor
    N.HD <- area/100   # Number of trees to consider for HD
    HT.HD <- plotdata$HT[order(plotdata$DBH, decreasing = TRUE)] # Sorting by diameters
    FT <- rep(0,length(plotdata$DBH))
    FT[1:ceiling(N.HD)] <- CF   # Note that this uses FT obtained from area.
    FT[ceiling(N.HD)] <- 100 - (sum(FT)-CF)  #The remaining proportion
    # #HD <- sum(HT.HD*FT)/100   #DBH times CF
    HD <- sum(HT.HD[1:ceiling(N.HD)]*FT[1:ceiling(N.HD)])/100   #DBH times CF
  }

  # Dealing with new columns

  if (sum(!is.na(extra))>0) {

    p <- length(extra)
    temp <- matrix(0,ncol=p, nrow=5)
    for (i in 1:p) {
      tr <- plotdata$FT*plotdata[,extra[i]]
      #colnames(tr)<-c('tr')
      #tr.sum <- data.frame(SP=plotdata$SPECIES, tr=tr)
      temp[1,i] <- sum(tr[plotdata$SPECIES==1], na.rm=TRUE)
      temp[2,i] <- sum(tr[plotdata$SPECIES==2], na.rm=TRUE)
      temp[3,i] <- sum(tr[plotdata$SPECIES==3], na.rm=TRUE)
      temp[4,i] <- sum(tr[plotdata$SPECIES==4], na.rm=TRUE)
      temp[5,i] <- sum(temp[,i])
    }
    colnames(temp)<-colnames(plotdata[extra])
    temp
    sd <- cbind(sd,temp)

  }

  return(list(sd=sd, DOM.SP=DOM.SP, PBAN=PBAN, PNHAN=PNHAN, BAind=plotdata$BAind, HD=HD))
}

# Note: need to define if this HD from tree data is used to replace the HD from site-dominant curves.
