% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site.R
\name{get_site}
\alias{get_site}
\title{Calculates missing stand-level variable from the set: AD, HD, SI}
\usage{
get_site(DOM.SP, ZONE, AD = NA, HD = NA, SI = NA)
}
\arguments{
\item{DOM.SP}{Dominant species (1: Rauli, 2: Roble, 3: Coigue)}

\item{ZONE}{Growth zone (1, 2, 3, 4)}

\item{AD}{Dominant age (year) of the stand}

\item{HD}{Dominant height (m) of dominant specie in the current stand}

\item{SI}{Site index (m) at reference dominant age of 20 years of the stand}
}
\value{
The missing stand-level parameter: AD, HD or SI
}
\description{
When two of stand level variables: dominat age (AD),
dominant height (HD) and site index (SI) are given, it returns the value of
the remaining stand-level parameter.
}
\examples{
# Example 1: Obtain Dominant Age
(AD <- get_site(DOM.SP=1, ZONE=2, HD=14, SI=10))

# Example 2: Obtain Dominant Height
(HD <- get_site(DOM.SP=1, ZONE=2, AD=25, SI=10))

# Example 3: Obtain Site Index
(SI <- get_site(DOM.SP=1, ZONE=2, AD=19, HD=13.5))
}
\references{
Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065, Chile

Gezan, S.A. and Moreno, P. (2000). Curvas de Sitio - Altura dominante para renovales
de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065, Chile
}
\author{
S. Gezan, S. Palmas and P. Moreno
}
