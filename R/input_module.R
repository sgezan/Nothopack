#' Module that process initial input information for future simulations/calculations
#'
#' It process input information for all future stand- or tree-level current
#' calculations or simulations from the user. It completes missing information that is required later.
#' It also makes the inventory processing of a given tree plot data.
#'
#' @param ZONE Growth zone (1, 2, 3, 4)
#' @param AD Dominant age (years) of the stand at present
#' @param AF Final dominant age (years) for simulation (if AD = AF no simulation is performed)
#' @param HD Dominant height (m) of dominant specie in the current stand
#' @param SI Site index (in m) at reference dominant age of 20 years of the stand
#' @param N Vector of number of trees (trees/ha) of the stand (1: Rauli, 2: Roble, 3: Coigue, 4:Others)
#' @param BA Vector of basal area (m2/ha) of the stand (1: Rauli, 2: Roble, 3: Coigue, 4:Others)
#' @param area Area of plot from tree-list data (m2) (required only for type = 'tree' or 'comp')
#' @param type Type of simulation required ('stand': stand-level, 'tree': tree-level,
#'           'comp': compatibility) (default = 'stand')
#' @param comptype Compatibility algorithm type ('PY': Proportional Yield, 'PG': Proportional Growth)
#'           (default = 'PY')
#' @param ddiam If TRUE a diameter distribution is generated (required only for type = 'stand')
#'           (default = FALSE)
#' @param NHA_model Number of model for N estimation (1:Original Reineke, 2:Re-fitted Reineke,
#'           3:Reineke with correction factor) (default = 3)
#' @param V_model Number of stand-level volume to use (1:with PNHAN, 2:without PNHAN)
#'           (required only for type = 'stand') (default = 1)
#' @param T_model Type of selected taper model (1:zone specific, 2:all zones) (default = 2)
#' @param IADBH_model Model to use for annual DBH increment (1, 2, 3, 4)
#'           (required only for type = 'tree' or 'comp') (default = 1)
#' @param ATHIN Age for thinning (years) (must be between AD and AF)
#' @param BARp Percentage of total basal area to remove (0-100) (required for ATHIN)
#' @param QDba Ratio of quadratic diameter of stand before against after thinning
#' @param FT.thin Vector of thinning decision for each tree belonging to the tree.list
#'           (0:to thin, 1:not to thin)
#' @param tree.list Tree-list for a plot with columns: ID, SPECIE, DBH, HT, SS, FT
#'           (required only for type = 'tree')
#'
#' @return List of data input and parameters to be traspassed for downstream modules (mainly core_module)
#' \itemize{
#'   \item \emph{Parameters} ZONE, DOM.SP, AD, HD, SI, SDIP, PBAN, PNHAN, AF, area, type, ddiam, comptype,
#'                NHA_model, V_model, T_model, IADBH_model, ATHIN, BARp, QD_ba,
#'                sp.table, tree.list, data.sim, DDist
#'   \item \emph{sp.table} Table with stand-level summary by SPECIES (1,2,3,4) and total (0)
#'                for N, BA, QD and VTHA
#'   \item \emph{tree.list} Data frame with tree list with missing values completed/estimated
#'                               (e.g. HT, FT, SS, VOL)
#'   \item \emph{data.sim} Data frame with summary information of yearly simulations
#'                               (including thinning)
#'    \item \emph{DDist} Data frame (matrix) with calculated (tree) or generated (stand)
#'                               diameter distribution
#'  }
#'
#' @author
#' S.A Gezan, S. Palmas and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example 1: Input stand-level data to predict volume all from stand-level functions
#' BA <- c(36.5, 2.8, 1.6, 2.4)
#' N <- c(464, 23, 16, 48)
#' input <- input_module(ZONE=1, AD=28, HD=18.5, AF=40, N=N, BA=BA, type='stand',
#'                       ddiam=FALSE, V_model=1)
#' input$SI
#' input$sp.table
#' input$data.sim
#'
#' # Example 2: Input stand-level data to predict volume by generating a diameter distribution
#' BA <- c(36.5, 2.8, 1.6, 2.4)
#' N <- c(464, 23, 16, 48)
#' stand <- input_module(ZONE=1, AD=28, HD=18.5, AF=28, N=N, BA=BA, type='stand',
#'                       ddiam=TRUE, T_model=2)
#' stand$sp.table
#' stand$data.sim
#' stand$DDist[5,,]
#' barplot(stand$DDist[5,,5], main='Diameter Distribution all species', xlab='DBH Class',
#'         names.arg=stand$DDist[5,,3], col='blue')
#'
#' # Example 3: Input plot-level data to process inventory (no simulation)
#' tree.data <- plot_example
#' head(tree.data,10)  # Original data
#' tree.inv <- input_module(ZONE=1, AD=28, type='tree', area=500, tree.list=tree.data, T_model=2)
#' tree.inv$DOM.SP    # Dominant specie
#' tree.inv$SI        # Site index
#' tree.inv$SDIP      # Stand density index
#' tree.inv$data.sim  # Summary data
#' tree.inv$sp.table  # Stand-level table by species
#' head(tree.inv$tree.list,10) # Completed tree-level data
#' plot(tree.inv$tree.list$DBH, tree.inv$tree.list$HT, ylab='HT (m)', xlab='DBH (cm)',
#'      xlim=c(0,40), ylim=c(0, 25), type='p', col='blue')
#' tree.inv$DDist[5,,] # Generated diameter distribution for all trees

input_module <- function(ZONE=NA,
                         AD=NA, AF=NA, HD=NA, SI=NA, N=NA, BA=NA, area=NA,
                         type='stand', comptype='PY', ddiam=FALSE,
                         NHA_model=3, V_model=1, T_model=2, IADBH_model=1,
                         ATHIN=NA, BARp=NA, QD_ba=NA, FT.thin=NA, tree.list=NA){

  DOM.SP <- NA
  SDIP <- NA
  PBAN <- NA
  PNHAN <- NA
  sdmatrix <- NA
  plotdata <- NA
  DDist <- NA

  ## Some traps for simulations (all types)

  # Errors with Age (AD, AF, ATHIN)
  if(is.na(AF)){
    warning('There is no final age for simulation ', call.=FALSE)
    AF<-AD
  }
  if(AD > AF){stop('The Final age of simulation should be larger (or equal) than the Initial Age')}

  # Checking details for compatibility
  if (type=='comp'){
      if(is.na(comptype)) {stop('Compatibility type (comptype) not indicated (PY or PG)')}
  }

  # Processing stand-level information
  if (type=='stand'){

    # Need for at least two of AD, HD, SI
    if(is.na(AD)==T && is.na(HD)==T && is.na(SI)==T |
       is.na(AD)==F && is.na(HD)==T && is.na(SI)==T |
       is.na(AD)==T && is.na(HD)==F && is.na(SI)==T |
       is.na(AD)==T && is.na(HD)==T && is.na(SI)==F){
      stop("Please provide information for AD, HD or SI (at least two of these).")
    }

    if(is.na(N)==T && is.na(BA)==T){
      stop('There is no enough information for simulations, provide with N and BA vectors')
    }

    N[5]<-sum(N)
    BA[5]<-sum(BA)
    QD<-N*0

    QD[1]<-get_stand(N=N[1], BA=BA[1])
    QD[2]<-get_stand(N=N[2], BA=BA[2])
    QD[3]<-get_stand(N=N[3], BA=BA[3])
    QD[4]<-get_stand(N=N[4], BA=BA[4])
    QD[5]<-get_stand(N=N[5], BA=BA[5])

    # Basic checking for dominant sp
    DOM.SP<-get_domsp(BA=BA[1:4])
    if (DOM.SP==99) {
      stop("This stand is not dominated by Nothofagus")
    }

    # Completing AD, HD or SI
    if (is.na(AD)){
      AD<-get_site(DOM.SP=DOM.SP, ZONE=ZONE, HD=HD, SI=SI)
    }
    if (is.na(HD)){
      HD<-get_site(DOM.SP=DOM.SP, ZONE=ZONE, AD=AD, SI=SI)
    }
    if (is.na(SI)){
      SI<-get_site(DOM.SP=DOM.SP, ZONE=ZONE, AD=AD, HD=HD)
      #print(SI)
    }

    # Completing the data frame for stand-table
    v1 <- c((1:4),0)
    v2 <- round(N,6)
    v3 <- round(BA,6)
    v4 <- round(QD,6)
    sdmatrix <- data.frame(cbind(v1,v2,v3,v4))
    names(sdmatrix) <- c('SPECIES','N','BA','QD')

    PBAN <- sum(BA[1:3])/(BA[5])
    PNHAN <- sum(N[1:3])/(N[5])

    # Calculation of SDI_percentage
    b1 <- 1.4112
    SDI <- N[5]*(QD[5]/25.4)^b1  # Good for all DOM.SP
    SDIMax <- NA
    if (DOM.SP==1 | DOM.SP==4) {
      SDIMax <- 1155.0 # Rauli and Mixed
    }
    if (DOM.SP==2) {
      SDIMax <- 908.8 # Roble
    }
    if (DOM.SP==3) {
      SDIMax <- 1336.9   # Coigue
    }
    SDIP <- round(100*SDI/SDIMax,6)

    # Calculates Stand-level volume (by stand-level equations)
    if (ddiam==FALSE) {

      # Calculating Stand-level volume (total)
      sdmatrix<-cbind(sdmatrix,VTHA=c(1:5)*0)
      if (V_model==1){
        VTHA <- Vmodule(BA=BA[5], HD=HD, PNHAN=PNHAN)
      } else {
        VTHA <- Vmodule(BA=BA[5], HD=HD)
      }
      sdmatrix[5,5]<-round(VTHA,6)

      # Assigning volume proportional to PBA
      PBA<-sdmatrix[1:4,3]/sdmatrix[5,3]
      sdmatrix[1:4,5]<-round(VTHA*PBA,6)

    }

    # Calculates Stand-level volume for simulations by diam.class (and generates diameter distributions)
    if (ddiam==TRUE) {

      # Generate stand-table from diameter distibution (does not contain volumes)
      stand.table<-diam_dist(sp.table=sdmatrix, HD=HD, DOM.SP=DOM.SP, ZONE=ZONE)
      m <- length(stand.table[1,,1])

      # Calculating class-level volume (total) by specie and class
      r_names<-c('DBH_ll','DBH_ul','D_class','H_class','N','BA','VT')
      DDist<-array(data=NA, dim=c(5,m,7), dimnames=list(c(1:5),c(1:m),r_names))
      DDist[,,-7]<-stand.table

      minN <- 0.1 # minimum Number of trees per diameter class to consider
      for (i in 1:m) {
        if (stand.table[1,i,5] < minN){ Vi.sp1<-0
        } else { Vi.sp1<-Vmodule_individual(SPECIES=1, ZONE=ZONE, DBH=stand.table[1,i,3],
                                            HT=stand.table[1,i,4], blength=stand.table[1,i,4], model=T_model) }
        if (stand.table[2,i,5] < minN){ Vi.sp2<-0
        } else { Vi.sp2<-Vmodule_individual(SPECIES=2, ZONE=ZONE, DBH=stand.table[2,i,3],
                                            HT=stand.table[2,i,4], blength=stand.table[2,i,4], model=T_model) }
        if (stand.table[3,i,5] < minN){ Vi.sp3<-0
        } else { Vi.sp3<-Vmodule_individual(SPECIES=3, ZONE=ZONE, DBH=stand.table[3,i,3],
                                            HT=stand.table[3,i,4], blength=stand.table[3,i,4], model=T_model) }
        if (stand.table[4,i,5] < minN){ Vi.sp4<-0
        } else { Vi.sp4<-Vmodule_individual(SPECIES=DOM.SP, ZONE=ZONE, DBH=stand.table[4,i,3],
                                            HT=stand.table[4,i,4], blength=stand.table[4,i,4], model=T_model)}
        DDist[1,i,7]<-round(Vi.sp1*stand.table[1,i,5],5)
        DDist[2,i,7]<-round(Vi.sp2*stand.table[2,i,5],5)
        DDist[3,i,7]<-round(Vi.sp3*stand.table[3,i,5],5)
        DDist[4,i,7]<-round(Vi.sp4*stand.table[4,i,5],5)
        DDist[5,i,7]<-DDist[1,i,7]+DDist[2,i,7]+DDist[3,i,7]+DDist[4,i,7]
      }

      # Assigning volume from generated stand-table
      sdmatrix<-cbind(sdmatrix,VTHA=c(1:5)*0)
      sdmatrix[1,5]<-round(sum(DDist[1,,7]),5)
      sdmatrix[2,5]<-round(sum(DDist[2,,7]),5)
      sdmatrix[3,5]<-round(sum(DDist[3,,7]),5)
      sdmatrix[4,5]<-round(sum(DDist[4,,7]),5)
      sdmatrix[5,5]<-round(sum(DDist[5,,7]),5)

    }

  }


  # Processing tree-level information
  if ((type=='tree') || (type=='comp')){

    # If area is missing
    if (is.na(area)) {stop('Area of plot not indicated')}

    #Check if the range of DBH is appropriate
    if (any(tree.list$DBH <= 5) | any(tree.list$DBH > 90)){
      warning('Some trees are smaller than 5 cm or larger than 90 cm DBH', call.=FALSE)
    }

    # No tree data from input
    if(nrow(tree.list)==0){
      stop('There is no tree data for type=tree')
    }

    # Completing FT if missing
    if (sum(is.na(tree.list$FT))>0) {
       FT<-matrix(1,nrow=nrow(tree.list),ncol=1)
       tree.list$FT<-FT*(10000/area)
    }

    # Completing HT if some are missing by regression
    if (sum(is.na(tree.list$HT))==length(tree.list$HT)) {stop('There is not data for tree heights, need to complete')}
    if (sum(!is.na(tree.list$HT))<length(tree.list$HT)) {
      warning('Some tree heights were estimated using model fitted with data (across all species)',call.=FALSE)
      warning('Verify that these are good estimates of total height.',call.=FALSE)
      tree.list$HT<-tree.ht(DBH=tree.list$DBH, HT=tree.list$HT, method=2)$HTFIN
    }

    # Adding DBH0 if missing
    if (is.null(tree.list$DBH0)) {
      tree.list$DBH0<-tree.list$DBH
    }

    # Completing SS if missing
    if (sum(is.na(tree.list$SS))>0){
      warning('Sociological Status (SS) was completed for all (or some) trees in plot if missing.',call.=FALSE)
      tree.list$SS<-tree_covariates(ID=tree.list$ID,FT=tree.list$FT,SPECIES=tree.list$SPECIES,DBH=tree.list$DBH,ZONE=ZONE)$SS
    }
    #if (sum(is.na(tree.list$SS))<length(tree.list$SS)) {stop('There is are some SS missing, need to complete')}

    # Collecting final summarized stand parameters (for sdmatrix)
    plotdata <- data.frame(tree.list$ID, tree.list$SPECIES,
                         tree.list$DBH, tree.list$HT,
                         tree.list$SS, tree.list$FT, tree.list$DBH0)
    colnames(plotdata) <- c('ID','SPECIES','DBH','HT','SS','FT','DBH0')
    params <- stand_parameters(plotdata=plotdata,area=area)
    sdmatrix <- params$sd
    PBAN <- params$PBAN
    PNHAN <- params$PNHAN

    # Using HD from tree data if is missing (with warning)
    if(is.na(HD)==T){
      HD<-params$HD
      warning('Dominant height (HD) missing, using the one obtained from the HT tree data', call.=FALSE)
      #print(HD)
    }

    # Basic checking for dominant sp
    DOM.SP<-params$DOM.SP
    if (DOM.SP==99) {
      warning("This stand is not dominated by Nothofagus", call.=FALSE)
      warning("set to be dominated by Nothofagus with greatest Basal area", call.=FALSE)
      DOM.SP<-which.max(sdmatrix$BA[1:3])
    }

    # Need for at least two of AD, HD, SI
    if(is.na(AD)==T && is.na(HD)==T && is.na(SI)==T |
       is.na(AD)==F && is.na(HD)==T && is.na(SI)==T |
       is.na(AD)==T && is.na(HD)==F && is.na(SI)==T |
       is.na(AD)==T && is.na(HD)==T && is.na(SI)==F){
      stop("Please provide information for AD, HD or SI (at least two of these).")
    }

    # Completing AD, HD or SI
    if (is.na(AD)){
      AD<-get_site(DOM.SP=DOM.SP, ZONE=ZONE, HD=HD, SI=SI)
    }
    if (is.na(HD)){
      HD<-get_site(DOM.SP=DOM.SP, ZONE=ZONE, AD=AD, SI=SI)
    }
    if (is.na(SI)){
      SI<-get_site(DOM.SP=DOM.SP, ZONE=ZONE, AD=AD, HD=HD)
      #print(SI)
    }

    # Collecting final stand parameters
    sdmatrix <- params$sd

    PBAN <- sum(sdmatrix[1:3,3])/sdmatrix[5,3]
    PNHAN <- sum(sdmatrix[1:3,2])/sdmatrix[5,2]

    # Calculation of SDI_percentage
    b1 <- 1.4112
    SDI <- sdmatrix[5,2]*(sdmatrix[5,4]/25.4)^b1  # Good for all DOM.SP
    if (DOM.SP==1 | DOM.SP==4) {
      SDIMax <- 1155.0 # Rauli and Mixed
    }
    if (DOM.SP==2) {
      SDIMax <- 908.8 # Roble
    }
    if (DOM.SP==3) {
      SDIMax <- 1336.9   # Coigue
    }
    SDIP <- round(100*SDI/SDIMax,6)

    ############
    # Calculates Stand-level volume (and obtains the diam. distribution)
    DDist<-NA
    m <- nrow(tree.list)

    # Estimating individual volume for each tree - place to add T_model
    # Also to estimate other product volumes (changing blength)
    tree.list$VIND <- mapply(Vmodule_individual, SPECIES=tree.list$SPECIES,
                                 ZONE=ZONE, DBH=tree.list$DBH,
                                 HT=tree.list$HT, blength=tree.list$HT, model=T_model)

    # Definition of diameter groups
    class <- 5
    diam <- seq(from=5,to=90,by=class)   # Diameter classes
    breaks. <- seq(from=7.5, to=87.5, by=5)  # improved Dclass calculation. Faster
    Dclass <- as.numeric(cut(tree.list$DBH,diam, right=FALSE))

    for (i in 1:m) {
      tree.list$Dclass[i] <- breaks.[Dclass[i]]
    }
    #View(tree.list)

    # Assigning trees to each diameter class (Dclass)
    #tree.list$Dclass <- breaks.[findInterval(x = tree.list$DBH,vec = breaks.,all.inside = TRUE)]
    tree.list$Nst <- tree.list$FT  #is it needed? Duplicate?
    tree.list$BAst <- tree.list$FT*(pi/4)*(tree.list$DBH^2/10000)
    tree.list$VOLst <- tree.list$FT*tree.list$VIND

    HT.agg <-aggregate(HT~Dclass+SPECIES,FUN=mean,data=tree.list)
    N.agg  <-aggregate(Nst~Dclass+SPECIES,FUN=sum,data=tree.list)
    BA.agg <-aggregate(BAst~Dclass+SPECIES,FUN=sum,data=tree.list)
    VT.agg <-aggregate(VOLst~Dclass+SPECIES,FUN=sum,data=tree.list)
    data.agg<-data.frame(HT.agg,N.agg[3],BA.agg[3],VT.agg[3])

    nclass <- length(diam)-1
    DBH_LL <- matrix(data=0,nrow=nclass,ncol=1)
    DBH_UL <- matrix(data=0,nrow=nclass,ncol=1)
    Dclass <- matrix(data=0,nrow=nclass,ncol=1)
    BAclass <- matrix(data=0,nrow=nclass,ncol=1)
    Hclass <- matrix(data=0,nrow=nclass,ncol=1)
    Vclass <- matrix(data=0,nrow=nclass,ncol=1)

    for (j in 1:(nclass)){
      DBH_LL[j] <- diam[j]                 # cm
      DBH_UL[j] <- diam[j+1]               # cm
      Dclass[j] <- (diam[j]+diam[j+1])/2   # cm
    }
    head<-data.frame(DBH_LL,DBH_UL,Dclass)

    r_names<-c('DBH_ll','DBH_ul','D_class','H_class','N','BA','VT')
    DDist<-array(data=0, dim=c(5,(length(diam)-1),7),
                 dimnames=list(c(1:5),c(1:(length(diam)-1)),r_names))

    data.agg1<-data.agg[which(data.agg$SPECIES==1),]
    Dclass1<-merge(head,data.agg1,by="Dclass",all=TRUE)
    data.agg2<-data.agg[which(data.agg$SPECIES==2),]
    Dclass2<-merge(head,data.agg2,by="Dclass",all=TRUE)
    data.agg3<-data.agg[which(data.agg$SPECIES==3),]
    Dclass3<-merge(head,data.agg3,by="Dclass",all=TRUE)
    data.agg4<-data.agg[which(data.agg$SPECIES==4),]
    Dclass4<-merge(head,data.agg4,by="Dclass",all=TRUE)

    Dclass1<-data.frame(head,HT=round(Dclass1[,5],2),N=round(Dclass1[,6],2),BA=round(Dclass1[,7],2),VT=round(Dclass1[,8],3))
    Dclass2<-data.frame(head,HT=round(Dclass2[,5],2),N=round(Dclass2[,6],2),BA=round(Dclass2[,7],2),VT=round(Dclass2[,8],3))
    Dclass3<-data.frame(head,HT=round(Dclass3[,5],2),N=round(Dclass3[,6],2),BA=round(Dclass3[,7],2),VT=round(Dclass3[,8],3))
    Dclass4<-data.frame(head,HT=round(Dclass4[,5],2),N=round(Dclass4[,6],2),BA=round(Dclass4[,7],2),VT=round(Dclass4[,8],3))
    Dclass5<-Dclass4

    for (j in 1:(nclass)){
      if (is.na(Dclass1$HT[j])) {
        Dclass1$HT[j]=0
        Dclass1$N[j]=0
        Dclass1$BA[j]=0
        Dclass1$VT[j]=0
      }
      if (is.na(Dclass2$HT[j])) {
        Dclass2$HT[j]=0
        Dclass2$N[j]=0
        Dclass2$BA[j]=0
        Dclass2$VT[j]=0
      }
      if (is.na(Dclass3$HT[j])) {
        Dclass3$HT[j]=0
        Dclass3$N[j]=0
        Dclass3$BA[j]=0
        Dclass3$VT[j]=0
      }
      if (is.na(Dclass4$HT[j])) {
        Dclass4$HT[j]=0
        Dclass4$N[j]=0
        Dclass4$BA[j]=0
        Dclass4$VT[j]=0
      }
    }

    Dclass5$N=Dclass1$N+Dclass2$N+Dclass3$N+Dclass4$N
    Dclass5$BA=Dclass1$BA+Dclass2$BA+Dclass3$BA+Dclass4$BA
    Dclass5$VT=Dclass1$VT+Dclass2$VT+Dclass3$VT+Dclass4$VT
    Dclass5$HT=round((Dclass1$N*Dclass1$HT+Dclass2$N*Dclass2$HT+Dclass3$N*Dclass3$HT+Dclass4$N*Dclass4$HT)/Dclass5$N,2)
    Dclass5$HT[is.na(Dclass5$HT)]<-0

    DDist[1,,]<-as.matrix(Dclass1)  #1: Rauli
    DDist[2,,]<-as.matrix(Dclass2)  #2: Roble
    DDist[3,,]<-as.matrix(Dclass3)  #3: Coigue
    DDist[4,,]<-as.matrix(Dclass4)  #4: Others
    DDist[5,,]<-as.matrix(Dclass5)  #0: Total

    # Assigning volume from generated stand-table
    sdmatrix<-cbind(sdmatrix,VTHA=c(1:5)*0)
    sdmatrix[1,5]<-round(sum(DDist[1,,7]),5)
    sdmatrix[2,5]<-round(sum(DDist[2,,7]),5)
    sdmatrix[3,5]<-round(sum(DDist[3,,7]),5)
    sdmatrix[4,5]<-round(sum(DDist[4,,7]),5)
    sdmatrix[5,5]<-round(sum(DDist[5,,7]),5)

    tree.list<-tree.list[,c(1:8)]

  }
#QD_ba=1
  data.sim<-data.frame(AGE=AD,HD=HD,NHA=sdmatrix[5,2],QD=sdmatrix[5,4],BA=sdmatrix[5,3],
                       NHAN=sum(sdmatrix[1:3,2]),NHA99=sdmatrix[4,2],BAN=sum(sdmatrix[1:3,3]),
                       BA99=sdmatrix[4,3],PBAN=PBAN,PNHAN=PNHAN,SI=SI,SDIP=SDIP,
                       VTOT=sdmatrix[5,5])
  #print(data.sim)
  # List that is output from here input somewhere else
  output <- list(ZONE=ZONE, DOM.SP=DOM.SP, AD=AD, HD=HD, SI=SI, SDIP=SDIP, PBAN=PBAN, PNHAN=PNHAN, AF=AF,
                area=area, type=type, ddiam=ddiam, comptype=comptype,
                NHA_model=NHA_model, V_model=V_model, T_model=T_model, IADBH_model=IADBH_model,
                ATHIN=ATHIN, BARp=BARp, QD_ba=QD_ba, FT.thin=FT.thin,
                sp.table=sdmatrix, tree.list=tree.list, data.sim=data.sim, DDist=DDist)

  return(input=output)

}
