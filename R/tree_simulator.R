#' Simulates individual-tree growth in diameter (DBH) and mortality for next year
#'
#' It simulates tree-level growth, mortality (and recruitment) of a given
#' stand requiring stand-level and tree-level parameters. Simulations are done using tree-level
#' models starting from intial age (AD0) until final age (ADF) in increments of 1 year.
#' The model of mortality is the same as whole-stand simulator.
#' Note, recruitment is currently not implemented.
#'
#' @param core.tree An object completely originated from core_module
#'
#' @return A series of elements and parameters with updated tables. Outputs are the same as input_module
#'         but these are updated for simulations. The main table is sp.table and tree.list that
#'         contain a stand table for the new simulations, and the tree list for the latest year, all
#'         based on the given parameters. Also, all other information is updated for time AF.
#'
#' @author
#' S.A. Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Moreno, P.; Palmas, S.; Escobedo, F.; Cropper, W.; Gezan, S. Individual-tree diameter growth models
#' for mixed nothofagus second growth forests in southern chile. Forests 2017(8), 506
#'
#' @examples
#' # Example: Simulation from Input from tree-level data (or file)
#' tree.list <- plot_example # use available data in example
#' head(tree.list)
#' plot <- input_module(ZONE=1, AD=28, AF=29, type='tree', area=500, tree.list=tree.list, T_model=2)
#' sims <- tree_simulator(core.tree=plot)
#' head(sims$tree.list)
#' plot$sp.table
#' sims$sp.table
#' plot$data.sim
#' sims$data.sim

tree_simulator <- function(core.tree=NULL){

  HT <- core.tree$tree.list$HT
  ZONE <- core.tree$ZONE
  DOM.SP <- core.tree$DOM.SP
  SI <- core.tree$SI
  HD <- core.tree$HD
  AD <- core.tree$AD
  IADBH_model <- core.tree$IADBH_model

  FTv<-rep(1,length(core.tree$tree.list$ID))    # simple vector
  # Estimating stand and individual variables from tree data
  input.data1<-tree_covariates(ID=core.tree$tree.list$ID, FT=core.tree$tree.list$FT,
                               SPECIES=core.tree$tree.list$SPECIES, DBH=core.tree$tree.list$DBH,
                               ZONE=ZONE, SS=core.tree$tree.list$SS)

  # Initial variables. This are updated a single year
  Fap <- input.data1$FT
  QDp <- input.data1$QD
  DBHp <- input.data1$DBH
  BALcp <- input.data1$BALc
  Ap <- FTv*core.tree$AD
  SSp <- input.data1$SS
  DAp <- FTv*core.tree$AD
  SSCALp <- input.data1$PScal
  SP <- input.data1$SPECIES
  SDIp <- input.data1$SDI  # provided from tree_covariates

    # ANNUAL INCREMENT (in mm) (working as a vector)
    # Note that all trees have the same age
    Gest <- FTv*0
    for (i in 1:length(Gest)) {
      Gest[i] <- AIDBH_module(SP=SP[i], DBH=DBHp[i], A=Ap[i], SS=SSp[i], SSCAL=SSCALp[i], BALc=BALcp[i],
                         ZONE=ZONE, AD=AD, SDI=SDIp[1], model=IADBH_model)  # model=3 earlier
    }
    DBHp1 <- DBHp+(Gest/10)   # all in cm
    # Dif.DBH <- Gest/10
    DBH0 <- DBHp

    # MORTALITY
    N1 <- NHAmodule(NHA0=input.data1$NHA, QD0=QDp, DOM.SP=DOM.SP, model=core.tree$NHA_model)  # How many survive (done as vector)
    Mortality <- unique(input.data1$NHA)-N1    # How many should die
    BALr <- SSCALp/sum(SSCALp)                 # Relative BAL for each tree
    Mi <- BALr*Mortality                       # Applying mortality based on BAL (PSCALp)
    Fap1 <- Fap-Mi                             # Adjusted expansion factor. The sum is equal to N1.

    # Obtaining NEW covariates using adjusted DBHp1 and Fap1
    data.temp <- tree_covariates(ID=core.tree$tree.list$ID, FT=Fap1, SPECIES=SP, DBH=DBHp1, ZONE=ZONE, SS=SSp)

    #UPDATE VALUES
    DBHp <- data.temp$DBH
    Fap <- data.temp$FT
    BALcp <- data.temp$BALc
    SDIp <- data.temp$SDI
    Ap <- Ap+1    # new age
    DAp <- DAp+1  # new dominant age
    SSCALp <- data.temp$PScal
    NHAp <- data.temp$NHA
    BAp <- data.temp$BA
    QDp <- data.temp$QD
    SSp <- data.temp$SS

  plot.proy <- data.frame(data.temp,HT,Fap,Ap,DAp,DBH0)

  # Height Increment Module
  HT0 <- HT
  m <- length(HT0)
  HTparam.0 <- matrix(data=0,nrow=m,ncol=1)
  HTparam.n <- matrix(data=0,nrow=m,ncol=1)
  Dif.HT <- matrix(data=0,nrow=m,ncol=1)
  HDp <- get_site(DOM.SP=DOM.SP, ZONE=ZONE, AD=DAp[1], SI=SI)

  #(HT<-height_param(DOM.SP=2, ZONE=2, HD=15, QD=12, DBH=24))

  # predicted height at the start (AD)
  HTparam.0 <- height_param(DOM.SP=DOM.SP, ZONE=ZONE,
                            HD=HD, QD=input.data1$QD[1],
                            DBH=input.data1$DBH)

  # predicted height at the end (AF)
  HTparam.n <- height_param(DOM.SP=DOM.SP, ZONE=ZONE,
                            HD=HDp, QD=QDp[1],
                            DBH=DBHp)

  # Difference in height growth
  Dif.HT <- HTparam.n-HTparam.0
  #Dif.HT <- Dif.HT %>% replace(.<0, 0)
  Dif.HT <- replace(Dif.HT,Dif.HT<0,0)

  # Assigning new height
  HTn <- HT0+Dif.HT
  #print(data.frame(HT0,HTparam.0,HTparam.n,HTn,Dif.HT))

  # Generating new tree.list
  tree.list <- data.frame(plot.proy$ID,
                       plot.proy$SPECIES,
                       round(plot.proy$DBH,6),
                       round(HTn,6),
                       round(plot.proy$SS,6),
                       round(plot.proy$FT,6),
                       round(plot.proy$DBH0,6))

  colnames(tree.list)<-c('ID','SPECIES','DBH','HT','SS','FT','DBH0')

  updated.plot <- input_module(ZONE=core.tree$ZONE,
                     AD=AD+1, AF=core.tree$AF+1, HD=HDp, SI=core.tree$SI,
                     #N=sp.table$N[1:4], BA=sp.table$BA[1:4],
                     area=core.tree$area, type=core.tree$type, ddiam=core.tree$ddiam, comptype=core.tree$comptype,
                     NHA_model=core.tree$NHA_model, V_model=core.tree$V_model, T_model=core.tree$T_model,
                     IADBH_model=core.tree$IADBH_model, # Originally 3
                     ATHIN=core.tree$ATHIN, BARp=core.tree$BARp, QD_ba=core.tree$QD_ba, FT.thin=core.tree$FT.thin,
                     tree.list=tree.list)

  return(output=updated.plot)

}

# Note: it calls again input_model at the end.
# Note: fixed to use AIDBH_model=3, a bit restrictive but true
