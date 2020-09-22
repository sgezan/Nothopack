#' General compatibility module between tree-level and stand-level simulations
#'
#' General compatibility module that takes tree-level simulations and
#' makes a compatibility with the stand-level simulations. Two compatibility methods are
#' available PY: proportional compatibility, and PG: proportional growth compatibility.
#' The input object (tree or stand) needs to have the option requesting compatibility
#' (i.e, type='comp') together with the method: comptype ('PY' or 'PG') (default = 'PY')
#'
#' @param sim.tree a simulated tree object from tree_simulator
#' @param sim.stand a simulated stand object from stand_simulator
#'
#' @author
#' S. Palmas, S.A. Gezan and P. Moreno
#'
#' @return A list of elements and parameters with updated tables with compatibility. Outputs are the same as
#'         input_module object but these are updated for simulations.
#'
#' @examples
#' # Example 1: Proportional yield (PY)
#' # Tree information
#' tree.list <- plot_example
#' plot.tree <- input_module(ZONE=2, AD=28, HD=23.5, AF=45, type='tree', area=500,
#'                           tree.list=tree.list, ddiam=FALSE)
#' sim.tree <- core_module(input=plot.tree)
#' head(sim.tree$tree.list)
#' # Stand information from plot.tree$sp.table
#' BA <- c(1.086,38.915,0.0,0.313)
#' N <- c(60,780,0,80)
#' plot.stand <- input_module(ZONE=2, AD=28, HD=23.5, AF=45, N=N, BA=BA,
#'                            type='stand', ddiam=FALSE)
#' sim.stand <- core_module(input=plot.stand)
#' # Summary tree and stand by specie table
#' sim.tree$sp.table
#' sim.stand$sp.table
#' # Requesting Compatibility, for method 'PY'
#' sim.tree$comptype = 'PY'
#' sim.tree$type = 'comp'
#' sim.comp <- comp_module(sim.tree=sim.tree, sim.stand=sim.stand)
#' sim.comp$sp.table
#'
#' # Example 2: Proportional growth (PG)
#' sim.tree$comptype = 'PG'
#' sim.tree$type = 'comp'
#' sim.comp <- comp_module(sim.tree=sim.tree, sim.stand=sim.stand)
#' sim.comp$sp.table


comp_module <- function(sim.tree=NA, sim.stand=NA){

    #p1.SIM <- sim.tree$comp.list$prob.surv        # Survival probability  ??
    NHA1.SIM <- sim.stand$sp.table$N[5]   # From stand_simulation
    BA.SIM <- sim.stand$sp.table$BA[5]    # Simulated BA from stand_simulator
    DBH1.SIM <- sim.tree$tree.list$DBH    # Simulated new diameter from tree_simulator
    DBH0 <- sim.tree$tree.list$DBH0       # Original DBH from tree_simulatior
    FT.SIM <- sim.tree$tree.list$FT       # Expansion factor from tree_simulation

    # Default option PY for comptype
    if (is.na(sim.tree$comptype)) { sim.tree$comptype <- 'PY' }

    # Proportional compatibility - PY
    if (sim.tree$comptype == 'PY'){
      FT.COMP <- FT.SIM*(NHA1.SIM/sum(FT.SIM))   # Adjusting FT (mortality) (fixed, moved earlier)
      DBH1.SIM.COMP <- sqrt( ((DBH1.SIM)^2)*(BA.SIM/(pi/40000))/(sum(FT.COMP*(DBH1.SIM)^2)) )  # Adjusting diameter
    }
    # Individual compatibility
    else if (sim.tree$comptype == 'PG'){   # It is not iterated for FT (see paper)

      #Adjusted FT1
      FT.nlm<-function(m) {
        abs(NHA1.SIM-sum(FT.SIM^m))
      }
      m.est<-optimize(FT.nlm, c(0.6,1.4))$minimum
      FT.COMP <- FT.SIM^m.est

      #Adjusted DBH1
      Num <- BA.SIM/(pi/40000)
      Den <- sum(FT.COMP*(DBH0^2))
      DBH1.SIM.COMP <- sqrt( DBH0^2 + (DBH1.SIM^2-DBH0^2)*(Num-Den)/(sum(FT.COMP*(DBH1.SIM^2))-sum(FT.COMP*(DBH0^2))) )
      #m <- NHA1.SIM/sum(FT.SIM)/sum(FT.SIM/NHA1.SIM) #does not affect. The same as tree_simulator then
      #m <- log(NHA1.SIM)/sum(log(FT.SIM))

    }

    sim.tree$tree.list$DBH <- DBH1.SIM.COMP  # Updated DBH with compatibility
    sim.tree$tree.list$FT <- FT.COMP         # Updated FT (mortality with compatibility)
    #print(head(sim.tree$tree.list))

    sim.comp <- input_module(ZONE=sim.tree$ZONE, AD=sim.tree$AD+1, AF=sim.tree$AF, HD=sim.tree$HD, SI=sim.tree$SI,
                             area=sim.tree$area, type=sim.tree$type, comptype=sim.tree$comptype,
                             ddiam=sim.tree$ddiam, NHA_model=sim.tree$NHA_model, V_model=sim.tree$V_model,
                             T_model=sim.tree$T_model, IADBH_model=sim.tree$IADBH_model,
                             ATHIN=sim.tree$ATHIN, BARp=sim.tree$BARp,
                             tree.list=sim.tree$tree.list)

    return(sim.comp=sim.comp)
}
