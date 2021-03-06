% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_module.R
\name{core_module}
\alias{core_module}
\title{Module that coordinates the communication between the different updating and simulation modules}
\usage{
core_module(input = NULL)
}
\arguments{
\item{input}{List created by input_module, that is used to pass all information to run simulations}
}
\value{
A series of elements and parameters with updated tables. Outputs are the same as input_module
        but these are updated for simulations. The main table updated is data.sim that contains a summary
        of the simulations from AD to AF for same parameters. Also, all other information is updated for
        time AF.
}
\description{
It coordinates the communication between the different critical modules of
updating and simulation including the stand- or tree-level information. It works as the core between input,
calculations, simulations and output. Some characteristics are:
1) it requires initially all elements that originate from an input module that is assumed to be complete,
2) input/output can be tree-list or stand parameters,
3) runs the simulation according to input parameters.
}
\examples{
# Example 1: Simulation using stand-level data - Volume calculated from stand functions
# (no generation diam. distr.)
BA <- c(36.5, 2.8, 1.6, 2.4)
N <- c(464, 23, 16, 48)
input <- input_module(ZONE=1, AD=28, HD=18.5, AF=40, N=N, BA=BA,
                      type='stand', ddiam=FALSE, V_model=1)
input
input$sp.table
input$DDist
input$data.sim

# Example 2: Proccesing stand-level data as inventory calculation generating
#            a diameter distribution
BA <- c(36.5, 2.8, 1.6, 2.4)
N <- c(464, 23, 16, 48)
stand <- input_module(ZONE=1, AD=28, HD=18.5, AF=28, N=N, BA=BA,
                      type='stand', ddiam=TRUE, T_model=1)
stand$sp.table
stand$DDist[5,,]

# Example 3: Simulation from Input from tree-level data (or file)
tree.list <- plot_example # use available data in example
head(tree.list)
plot <- input_module(ZONE=1, AD=28, AF=33, type='tree', area=500,
                     tree.list=tree.list, T_model=2)
plot$DDist[5,,]
head(plot$tree.list)
sims <- core_module(input=plot)
head(sims$tree.list)
plot$sp.table
sims$sp.table
plot$data.sim
sims$data.sim
data.temp <- as.matrix(data.frame(Initial=plot$DDist[5,,5], Final=sims$DDist[5,,5]))
barplot(data.temp, main='Diameter Distribution all species', xlab='DBH Class',
        beside=TRUE, col=cbind(rep('blue',17), rep('darkblue',17)))
}
\references{
Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
}
\author{
S.A. Gezan, S. Palmas and P. Moreno
}
