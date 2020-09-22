#' Module that coordinates the communication between the different updating and simulation modules
#'
#' It coordinates the communication between the different critical modules of
#' updating and simulation including the stand- or tree-level information. It works as the core between input,
#' calculations, simulations and output. Some characteristics are:
#' 1) it requires initially all elements that originate from an input module that is assumed to be complete,
#' 2) input/output can be tree-list or stand parameters,
#' 3) runs the simulation according to input parameters.
#'
#' @param input List created by input_module, that is used to pass all information to run simulations
#'
#' @return A series of elements and parameters with updated tables. Outputs are the same as input_module
#'         but these are updated for simulations. The main table updated is data.sim that contains a summary
#'         of the simulations from AD to AF for same parameters. Also, all other information is updated for
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
#' # Example 1: Simulation using stand-level data - Volume calculated from stand functions
#' # (no generation diam. distr.)
#' BA <- c(36.5, 2.8, 1.6, 2.4)
#' N <- c(464, 23, 16, 48)
#' input <- input_module(ZONE=1, AD=28, HD=18.5, AF=40, N=N, BA=BA,
#'                       type='stand', ddiam=FALSE, V_model=1)
#' input
#' input$sp.table
#' input$DDist
#' input$data.sim
#'
#' # Example 2: Proccesing stand-level data as inventory calculation generating
#' #            a diameter distribution
#' BA <- c(36.5, 2.8, 1.6, 2.4)
#' N <- c(464, 23, 16, 48)
#' stand <- input_module(ZONE=1, AD=28, HD=18.5, AF=28, N=N, BA=BA,
#'                       type='stand', ddiam=TRUE, T_model=1)
#' stand$sp.table
#' stand$DDist[5,,]
#'
#' # Example 3: Simulation from Input from tree-level data (or file)
#' tree.list <- plot_example # use available data in example
#' head(tree.list)
#' plot <- input_module(ZONE=1, AD=28, AF=33, type='tree', area=500,
#'                      tree.list=tree.list, T_model=2)
#' plot$DDist[5,,]
#' head(plot$tree.list)
#' sims <- core_module(input=plot)
#' head(sims$tree.list)
#' plot$sp.table
#' sims$sp.table
#' plot$data.sim
#' sims$data.sim
#' data.temp <- as.matrix(data.frame(Initial=plot$DDist[5,,5], Final=sims$DDist[5,,5]))
#' barplot(data.temp, main='Diameter Distribution all species', xlab='DBH Class',
#'         beside=TRUE, col=cbind(rep('blue',17), rep('darkblue',17)))


core_module <- function(input=NULL){

  AD <- input$AD
  AF <- input$AF
  type <- input$type
  ddiam <- input$ddiam
  comptype <- input$comptype
  ATHIN <- input$ATHIN
  BARp <- input$BARp
  QD_ba <- input$QD_ba
  FT.thin <- input$FT.thin
  sp.table <- NA
  tree.list <- NA
  data.sim <- input$data.sim
  DDiam <- NA

  # Some small initial checks
  if (AD==AF) {
    warning('Initial and final age of simulations are identical (AD=AF)', call.=FALSE)
  }

  #########################
  # Start main annual cycle

  AGE<-AD
  for (k in AD:AF) {

    # STAND simulation (no comp.)
    if (type=='stand') {

      if (!is.na(ATHIN)){
        if(ATHIN < AD){
          warning('The thinning age is preliminar to the final age of simulation.
                  Thinning time changed to initial age (ATHIN = AD).', call.=FALSE)
          ATHIN<-AD
        }
        if(ATHIN > AF){
          warning('The thinning age is posterior to the final age of simulation.
                  Thinning time changed to final age (ATHIN = AF).', call.=FALSE)
          ATHIN<-AF
        }
        if(is.na(BARp)){
          warning('Proportion of basal area to remove from thinning not indicated.
                  Hence, no thinning is implemented.', call.=FALSE)
          ATHIN<-NA
        }
        if(is.na(QD_ba)) {
          warning('Ratio of QD before/after not specified set to 1.0', call.=FALSE)
          QD_ba<-1.0
        }
        if(input$ddiam==TRUE) {
          input$ddiam <- FALSE
          warning('Thinning is not allowed with generation of diameter distribution. Simulation is run without generating diameter distribution.',call.=FALSE)
        }
      }

      # If there is thinning
      if (!is.na(ATHIN)){
        if (AGE==ATHIN) { # Are we in thinning age? Yes
          input$AD <- AGE
          sim.stand <- thin_module(core.thin=input)
          data.sim <- rbind(data.sim,sim.stand$data.sim)
          # Recalculate things
          #sim.stand <- stand_simulator(core.stand=input)
          sim.stand$data.sim <- data.sim
          input <- sim.stand
          sims <- sim.stand
        }
        if (AGE<AF) {  # We simulate if we still have years
          input$AD <- AGE
          input$AF <- AGE+1
          sim.stand <- stand_simulator(core.stand=input)
          data.sim <- rbind(data.sim,sim.stand$data.sim)
          sim.stand$data.sim <- data.sim
          input <- sim.stand
          sims <- sim.stand
          AGE<-AGE+1
        }
      }

      # If there is no thinning then we just simulate
      if (is.na(ATHIN) & AGE<AF){
        input$AD <- AGE
        input$AF <- AGE+1
        sim.stand <- stand_simulator(core.stand=input)
        data.sim <- rbind(data.sim,sim.stand$data.sim)
        sim.stand$data.sim <- data.sim
        input <- sim.stand
        sims <- sim.stand
        AGE<-AGE+1
      } else {
        sims <- input
      }

      if (ddiam==FALSE) {
        sims$DDist <- NA
        sims$tree.list <- NA
      }

    }

    # TREE simulation (no comp.)
    if (type=='tree') {

      if (!is.na(ATHIN)){
        if(ATHIN < AD){
          warning('The thinning age is preliminar to the final age of simulation.
                  Thinning time changed to initial age (ATHIN = AD).', call.=FALSE)
          ATHIN<-AD
        }
        if(ATHIN > AF){
          warning('The thinning age is posterior to the final age of simulation.
                  Thinning time changed to final age (ATHIN = AF).', call.=FALSE)
          ATHIN<-AF
        }
      }

      # If there is thinning
      if (!is.na(ATHIN)){
        if (AGE==ATHIN) { # Are we in thinning age? Yes
          input$AD <- AGE
          sim.tree <- thin_module(core.thin=input)
          # This line is added
          #sim.tree <- tree_simulator(core.tree=sim.tree)
          data.sim <- rbind(data.sim,sim.tree$data.sim)
          sim.tree$data.sim <- data.sim
          input <- sim.tree
          sims <- sim.tree
          #AGE<-AGE+1
        }
        if (AGE<AF) {  # We simulate if we still have years
          input$AD <- AGE
          input$AF <- AGE+1
          sim.tree <- tree_simulator(core.tree=input)
          data.sim <- rbind(data.sim,sim.tree$data.sim)
          sim.tree$data.sim <- data.sim
          input <- sim.tree
          sims <- sim.tree
          AGE<-AGE+1
        }
      }

      # If there is not thinning then we just simulate
      if (is.na(ATHIN) & AGE<AF){
        input$AD <- AGE
        input$AF <- AGE+1
        sim.tree <- tree_simulator(core.tree=input)
        data.sim <- rbind(data.sim,sim.tree$data.sim)
        sim.tree$data.sim <- data.sim
        input <- sim.tree
        sims <- sim.tree
        AGE<-AGE+1
      } else {
        sims <- input
      }

    }

    # COMPATIBILITY simulation
    if (type=='comp'& AGE<AF) {

      if(!is.na(input$ATHIN)) {
        input$ATHIN <- NA
        warning('Thinning is not allowed with compatibility. Simulation is run ignoring thinning',call.=FALSE)
      }

      input$AD <- k
      input$AF <- k+1

      # Stand simulation
      input$type <- 'stand'
      sim.stand <- core_module(input=input)
      input$type <- 'tree'
      sim.tree <- core_module(input=input)

      # Compatibility
      sim.comp <- comp_module(sim.tree=sim.tree, sim.stand=sim.stand)
      data.sim <- rbind(data.sim,sim.comp$data.sim)
      sim.comp$type <- 'comp'
      input <- sim.comp
      sims <- sim.comp
      AGE<-AGE+1

    }
    #data.sim <- rbind(data.sim,sim.comp$data.sim)
    sims$data.sim<-data.sim
    sims$AD <- AD
    sims$AF <- AF

  }

  return(output=sims)

}
