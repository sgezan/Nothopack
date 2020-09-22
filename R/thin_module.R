#' Executes thinning for a plot according to given stand- or tree-level information
#'
#' It executes a thinning for a plot based on stand- or tree-level information.
#' The thinning goes across all species (Nothofagus and other species) using identical
#' criteria.
#'
#' @param core.thin input object generated from \emph{core_module}
#'
#' @return A series of elements and parameters with updated tables. Outputs are the same as \emph{input_module}
#'         but these are updated for simulations. The main table is \emph{sp.table} that contains a stand table
#'         for the new simulations for the given parameters. Also, all other information is updated for
#'         time AF.
#'
#' @author
#' S.A. Gezan
#'
#' @examples
#' # Example 1: Stand thinning (equal rate for every diameter class) at 30 years
#' #    from initial stand-level data with a residual basal area porcentual of 50%
#' BA <- c(1.09, 38.92, 0, 0.31)
#' N <- c(60, 780, 0, 80)
#' stand.input <- input_module(ZONE=2, AD=28, AF=70, ATHIN=30, HD=18.4, N=N, BA=BA,
#'                             type='stand', ddiam=FALSE, BARp=50, QD_ba=1)
#' sims.thin <- core_module(input=stand.input)
#' sims.thin$data.sim
#' plot(sims.thin$data.sim$AGE, sims.thin$data.sim$BA, ylim=c(0,60),
#'      xlab='Age (years)', type='o', col='blue', ylab='BA (m2/ha)')
#' plot(sims.thin$data.sim$AGE, sims.thin$data.sim$VTOT, ylim=c(0,600),
#'      xlab='Age (years)', type='o', col='blue', ylab='VTOT (m3/ha)')
#'
#' # Example 2: Simulation with tree thinning at 30 years starting from initial tree-level
#' #    data (or file). FT.thin is a random vector of thinning, 0=out, 1=in (to be added)
#' tree.list <- plot_example
#' head(tree.list)
#' FT.thin <- round(0.3+runif(length(tree.list$FT),0))
#' head(FT.thin, 20)
#' plot.input <- input_module(ZONE=1, AD=28, AF=40, type='tree', ATHIN=30, FT.thin=FT.thin,
#'                            area=500, tree.list=tree.list, T_model=2)
#' head(plot.input$tree.list, 10)
#' sim.thin <- core_module(input=plot.input)
#' head(sim.thin$tree.list, 10)
#' plot.input$sp.table
#' sim.thin$sp.table
#' sim.thin$data.sim
#' plot(sim.thin$data.sim$AGE, sim.thin$data.sim$BA, ylim=c(0,60),
#'      xlab='Age (years)', type='o', col='blue', ylab='BA (m2/ha)')
#' plot(sim.thin$data.sim$AGE, sim.thin$data.sim$VTOT, ylim=c(0,400),
#'      xlab='Age (years)', type='o', col='blue', ylab='VTOT (m3/ha)')

thin_module <- function(core.thin=NULL){

  # STAND thinning
  if (core.thin$type=='stand') {

    # Initial information
    PBAN<-core.thin$PBAN
    PNHAN<-core.thin$PNHAN
    area<-core.thin$area

    BARp<-core.thin$BARp   # Percentage of total basal area to remove (0-100%) BAR = 30
    QD_ba<-core.thin$QD_ba
    if(is.na(QD_ba)) {
      warning('Ratio of QD before/after not specified set to 1.0', call.=FALSE)
      QD_ba<-1.0
    }
    sp.table<-core.thin$sp.table
    tree.list<-core.thin$tree.list

    if (QD_ba<0.7) {
      stop('Thinning from below too strong (0.7 < QD_ba < 1.3)')
    }
    if (QD_ba>=1.3) {
      stop('Thinning from above too strong (0.7 < QD_ba < 1.3)')
    }

    sp.table$BA <- sp.table$BA*(100-BARp)/100
    sp.table$QD <- sp.table$QD/QD_ba

    sp.table$N[1]<-round(get_stand(BA=sp.table$BA[1], QD=sp.table$QD[1]),5)
    sp.table$N[2]<-round(get_stand(BA=sp.table$BA[2], QD=sp.table$QD[2]),5)
    sp.table$N[3]<-round(get_stand(BA=sp.table$BA[3], QD=sp.table$QD[3]),5)
    sp.table$N[4]<-round(get_stand(BA=sp.table$BA[4], QD=sp.table$QD[4]),5)
    sp.table$N[5]<-round(get_stand(BA=sp.table$BA[5], QD=sp.table$QD[5]),5)

    output <- input_module(ZONE=core.thin$ZONE,
                           AD=core.thin$AD, AF=core.thin$AF, HD=core.thin$HD, SI=core.thin$SI,
                           N=sp.table$N[1:4], BA=sp.table$BA[1:4],
                           area=core.thin$area, type=core.thin$type, ddiam=core.thin$ddiam, comptype=core.thin$comptype,
                           NHA_model=core.thin$NHA_model, V_model=core.thin$V_model, T_model=core.thin$T_model,
                           IADBH_model=core.thin$IADBH_model,
                           ATHIN=core.thin$ATHIN, BARp=core.thin$BARp, QD_ba=core.thin$QD_ba, FT.thin=core.thin$FT.thin,
                           tree.list=core.thin$tree.list)

  }

  # TREE thinning
  if (core.thin$type=='tree') {

    if (length(core.thin$FT.thin)!=length(core.thin$tree.list$FT)) {
      warning('Length of vector of thinning decision is not the same as the list of trees. Thinning is not ocurring.', call.=FALSE)
      core.thin$FT.thin <- core.thin$tree.list$FT*1
    }

    core.thin$tree.list$FT<-core.thin$tree.list$FT*core.thin$FT.thin
    core.thin$tree.list<-core.thin$tree.list[core.thin$tree.list$FT>0,]

    # This is the time to update sp.table
    output <- input_module(ZONE=core.thin$ZONE,
                         AD=core.thin$AD, AF=core.thin$AF, HD=core.thin$HD, SI=core.thin$SI,
                         #N=sp.table$N[1:4], BA=sp.table$BA[1:4],
                         area=core.thin$area, type=core.thin$type, ddiam=core.thin$ddiam, comptype=core.thin$comptype,
                         NHA_model=core.thin$NHA_model, V_model=core.thin$V_model, T_model=core.thin$T_model,
                         IADBH_model=core.thin$IADBH_model,
                         ATHIN=core.thin$ATHIN, BARp=core.thin$BARp, QD_ba=core.thin$QD_ba, FT.thin=core.thin$FT.thin,
                         tree.list=core.thin$tree.list)
    #print(output)

    QD_ba<-core.thin$sp.table[5,4]/output$sp.table[5,4]
    BARp<-100*(core.thin$sp.table[5,3]-output$sp.table[5,3])/core.thin$sp.table[5,3]
    print(paste0('Tree thinning resulted in QD_ba = ', round(QD_ba,3)))
    print(paste0('Tree thinning resulted in BARp = ', round(BARp,2)))
    if (QD_ba<0.7) {
      warning('Thinning from below too strong', call.=FALSE)
    }
    if (QD_ba<0.7) {
      warning('Thinning from below too strong', call.=FALSE)
    }
    if (QD_ba>=1.3) {
      warning('Thinning from above too strong', call.=FALSE)
    }

  }

  return(output)

}

# Note
# - DDiam is generated as if it is not thinned (all D_classes)
# - ddiam=TRUE and ATHIN can not occur (no generation of ddiam if there is thinning)
