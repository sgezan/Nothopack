% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/height_param.R
\name{height_param}
\alias{height_param}
\title{Estimates total individual tree height in function of its DBH and stand-level parameters}
\usage{
height_param(DOM.SP, ZONE, HD = NA, QD = NA, DBH = NA)
}
\arguments{
\item{DOM.SP}{Dominant species (1: Rauli, 2: Roble, 3: Coigue) of the stand}

\item{ZONE}{Growth zone (1, 2, 3, 4) of the stand}

\item{HD}{Dominant height (m) of dominant specie in the current stand}

\item{QD}{Quadratic diameter (cm) of the stand}

\item{DBH}{Diameter at breast height (cm) of tree}
}
\value{
Individual total tree height HT (m)
}
\description{
Estimates individual tree height form a model of stand-level parameters:
dominant height (HD, m), quadratic diameter (QD, cm), and diameter at breast height (DBH, cm).
}
\examples{
# Example: Estimates HT for trees with DBH=24 for different stands
(HT <- height_param(DOM.SP=2, ZONE=2, HD=15, QD=12, DBH=24))
(HT <- height_param(DOM.SP=2, ZONE=2, HD=15, QD=18, DBH=24))
(HT <- height_param(DOM.SP=2, ZONE=2, HD=19, QD=12, DBH=24))
(HT <- height_param(DOM.SP=2, ZONE=2, HD=19, QD=18, DBH=24))
}
\references{
Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
}
\author{
S. Palmas, S.A. Gezan and P. Moreno
}
