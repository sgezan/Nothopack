#' Obtains/calculates stand- and tree-level covariates that feed into the individual DBH tree growth model
#'
#' It obtains/calculates several stand- and tree-level covariables
#' that form part of the individual-tree DBH growth model (AIDBH). The input corresponds to tree data
#' originated from an inventory of a single plot, and they should all be on the same tree ID order.
#' Where there is missing data (e.g, SS) this information is completed.
#'
#' @param ID Vector of identification of each tree
#' @param FT Vector of expansion factor (\emph{i.e.}, FT = 1/area) for each tree
#' @param SPECIES Vector of species code (1:Rauli, 2:Roble, 3:Coigue, 4:Other) for each tree
#' @param DBH Vector of diameter breast height (DBH, cm) for each tree
#' @param ZONE Growth zone (1, 2, 3, 4) for the plot
#' @param SS Vector of sociological status (range 1-4, where 1:Dominant, 2:Codominant, 3:Intermediate, 4:Supressed)
#'
#' @return A data frame containing the following vectors:
#'          ID (tree id), FT (expansion factor), DBH (diameter breast height, cm),
#'          ba (basal area tree, m2/ha), bac (basal area Nothofagus cohort, m2/ha),
#'          SPECIES (tree specie), ZONE (growth zone), SPZONA (concatenation specie-zone),
#'          N (number of trees/ha), BA (basal area, m2/ha), QD (quadratic diameter, cm),
#'          SDI (stand density index, trees/ha), SS (sociological status, range 1-4),
#'          BAL (basal area of larger trees, m2/ha), BALc (basal area of larger trees cohort Nothofagus, m2/ha),
#'          PSCAL (calculated sociological status as BAL/BA)
#'
#' @author
#' P. Moreno and S.A. Gezan
#'
#' @references
#' Moreno et al. (2017). Individual-tree diameter growth models for mixed Nothofagus
#' second growth forests in southern Chile. Forests 8(12), 506.
#'
#' @examples
#' # Example 1: Simple data with SS to complete
#' ID <- c(1, 2, 3, 4)
#' FT <- c(400, 400, 400, 400)
#' SPECIES <- c(1, 1, 2, 4)
#' DBH <- c(15, 25, 37, 20)
#' SS <- c(NA, 2, NA, 2)
#' tree.fill <- tree_covariates(ID=ID, FT=FT, SPECIES=SPECIES, DBH=DBH, ZONE=1, SS=SS)
#' tree.fill
#'
#' # Example 2: Complete inventory full plot
#' tree.list <- plot_example
#' tree.list$FT <- 10000/500
#' covar <- tree_covariates(ID=tree.list$ID, FT=tree.list$FT, SPECIES=tree.list$SPECIES,
#'                          DBH=tree.list$DBH, ZONE=2, SS=tree.list$SS)
#' head(covar, 10)

tree_covariates<-function(ID=NA, FT=NA, SPECIES=NA, DBH=NA, ZONE=NA, SS=NA){

  N <- length(ID)                              # Getting the number of trees in the plot
  NHA <- sum(FT)                               # Getting the number of trees in the stand by hectare
  BA <- sum(pi*(DBH/2/100)^2*FT, na.rm = TRUE) # Total basal area by hectare
  QD <- (100*((4/pi)*(BA/NHA))^0.5)            # Quadratic Diameter,
  SDI<- NHA*((25.4/QD)^(-1.4112))              # Stand density index
  SPZONA<-SPECIES*10+as.numeric(ZONE)          # Why times a zone?

  # BAL & BALc calculation
  ba<-pi*((DBH)^2)/40000  # in m2
  bac<-0

  for (fila in (1:N)) {
    bac[fila] <- if (SPECIES[fila] != 4) {ba[fila]} else {0}
  }

   #Creating a table BAL and BALc info
   # Temp.data<-data.frame(ID,FT,DBH,ba,bac,SPECIES,ZONE,SPZONA,NHA,BA,QD,SDI,SS) %>%
   #   arrange(desc(ba)) %>%
   #   mutate(
   #     BAL = (lag(ba) %>% replace(is.na(.), 0) %>% cumsum) * FT,
   #     BALc = (lag(bac) %>% replace(is.na(.), 0) %>% cumsum) * FT,
   #     BALc = BALc %>% replace(is.na(.), mean(BALc)),   #replacing to the mean
   #     PScal = BAL/BA
   #   )

  # Creating a table with BAL, BALc and PScal
  Temp.data <- data.frame(ID,FT,DBH,ba,bac,SPECIES,ZONE,SPZONA,NHA,BA,QD,SDI,SS)
  Temp.data <- Temp.data[order(Temp.data$ba,decreasing=TRUE),]  # sorting by ba
  Temp.data$BAL <- cumsum(replace(lag(Temp.data$ba),is.na(lag(Temp.data$ba)),0))*FT
  Temp.data$BALc <- cumsum(replace(lag(Temp.data$bac),is.na(lag(Temp.data$bac)),0))*FT
  Temp.data$BALc <- replace(Temp.data$BALc,is.na(Temp.data$BALc),mean(Temp.data$BALc))
  Temp.data$PScal <- Temp.data$BAL/Temp.data$BA

  # Sociological status calculated by competence and stocking (completed where is missing)
  for (i in (1:nrow(Temp.data))){
    if (is.na(Temp.data$SS[i])){
      #Temp.data$SS[i]<-tree.SS(Temp.data$PScal[i])
      Temp.data$SS[i]<-1+Temp.data$PScal[i]*3
    }
  }

  # Needs to go back to original order of trees
  Temp.data<-Temp.data[order(Temp.data$ID),]

  return(Temp.data)

}
