#' Incorporates new young trees (ingrowth/recruitment) to a current stand's
#' diameter distribution (module in development)
#'
#' It adds a 1% of the total number of trees (N, trees/ha), which is
#' assigned to the species proportional to their number of trees. These are then
#' assigned to the first diameter class and reported to be added.
#' In the future it should estimate the 'space' for new trees (high competition, low number of
#' new trees). The procedure generates a new diameter distribution and then add trees to small classes.
#' NOTE 1: This model is not fully implemented nor incorporated into the simulator.
#' NOTE 2: Adjustment of VTHA is not included here.
#'
#' @param sp.table table with stand-level information by specie and total, columns are
#'            SPECIES, N, BA, QD, and table is grouped by SPECIES (1:Rauli, 2:Roble, 3:Coigue, 4:Other, 0:Total).
#' @param HD Dominant height (m) of dominant specie in the current stand
#' @param DOM.SP Dominant specie (1: Rauli, 2: Roble, 3: Coigue)
#' @param ZONE Growth zone (1, 2, 3, 4)
#'
#' @return A list with: Dd (data frame with updated diameter distribution including
#' recruits (N and BA are updated); Rec.SP (vector or random proportional assignation
#' of recruits to each specie; and sp.table (data frame with updated sp.table including
#' recuits.
#'
#' @author
#' S. Palmas, S.A. Gezan and P. Moreno
#'
#' @references
#' Gezan, S.A. and Ortega, A. (2001). Desarrollo de un Simulador de Rendimiento para
#' Renovales de Roble, Rauli y Coigue. Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example: Evaluating a new recruitment
#' BA <- c(36.5, 2.8, 1.6, 2.4)
#' N <- c(464, 23, 16, 48)
#' input <- input_module(ZONE=1, AD=28, HD=18.5, N=N, BA=BA, type='stand')
#' Dd <- diam_dist(sp.table=input$sp.table, HD=input$HD, DOM.SP=input$DOM.SP, ZONE=input$ZONE)
#' Dd.R <- recruit_module(sp.table=input$sp.table, HD=input$HD, DOM.SP=input$DOM.SP,
#'                        ZONE=input$ZONE)
#' Dd.R$Rec.SP
#' input$sp.table
#' Dd.R$sp.table
#' head(Dd[1,,])
#' head(Dd.R$Dd[1,,])

recruit_module <- function(sp.table=NA, HD=NA, DOM.SP=NA, ZONE=NA){

  vN <- sp.table$N[1:4]

  # Generates the diametric distribution of the current stand
  Dd <- diam_dist(sp.table=sp.table, HD=HD, DOM.SP=DOM.SP, ZONE=ZONE)

  # Estimates the probabilitt based on the current number of trees per species
  # Maybe it should be changed to some other function. Maybe based on dominant species, and density.
  # This is because some species may be recruited differently on different conditions
  # (e.g. rauli may stand shading better so it will be recruited better at higher densities)
  prob <- vN/sum(vN, na.rm = TRUE) # Probabilities of recruitment
  #print(prob)

  # This is the number of species to recruit
  #Nrecruit <- 3 # should be a function of BA, N and HD
  Nrecruit <- 0.01*sp.table$N[5]
  #print(Nrecruit)

  #recruited trees species
  Rec.SP <- table(factor(sample(x = c(1,2,3,4), size=Nrecruit, replace=TRUE, prob=prob),  # A list of n.trees species
                    levels = c(1,2,3,4)))

  # Updating the diametric distribution - N
  Dd[1,1,5] <- Dd[1,1,5] + Rec.SP[1]
  Dd[2,1,5] <- Dd[2,1,5] + Rec.SP[2]
  Dd[3,1,5] <- Dd[3,1,5] + Rec.SP[3]
  Dd[4,1,5] <- Dd[4,1,5] + Rec.SP[4]
  Dd[5,1,5] <- Dd[1,1,5] + Dd[2,1,5] + Dd[3,1,5] + Dd[4,1,5]

  # Updating the diametric distribution - BA
  Dd[1,1,6] <- round( Dd[1,1,5]*(pi/4)*(Dd[1,1,3]/100)^2, 5)
  Dd[2,1,6] <- round( Dd[2,1,5]*(pi/4)*(Dd[2,1,3]/100)^2, 5)
  Dd[3,1,6] <- round( Dd[3,1,5]*(pi/4)*(Dd[3,1,3]/100)^2, 5)
  Dd[4,1,6] <- round( Dd[4,1,5]*(pi/4)*(Dd[4,1,3]/100)^2, 5)
  Dd[5,1,6] <- round( Dd[5,1,5]*(pi/4)*(Dd[5,1,3]/100)^2, 5)

  # Updating sp.table
  update.sp.table <- sp.table
  update.sp.table[1,2] <- sum(Dd[1,,5])  #N
  update.sp.table[2,2] <- sum(Dd[2,,5])  #N
  update.sp.table[3,2] <- sum(Dd[3,,5])  #N
  update.sp.table[4,2] <- sum(Dd[4,,5])  #N
  update.sp.table[5,2] <- sum(Dd[5,,5])  #N
  update.sp.table[1,3] <- sum(Dd[1,,6])  #BA
  update.sp.table[2,3] <- sum(Dd[2,,6])  #BA
  update.sp.table[3,3] <- sum(Dd[3,,6])  #BA
  update.sp.table[4,3] <- sum(Dd[4,,6])  #BA
  update.sp.table[5,3] <- sum(Dd[5,,6])  #BA
  update.sp.table[1,4] <- get_stand(BA=update.sp.table[1,3], N=update.sp.table[1,2])  #QD
  update.sp.table[2,4] <- get_stand(BA=update.sp.table[2,3], N=update.sp.table[2,2])  #QD
  update.sp.table[3,4] <- get_stand(BA=update.sp.table[3,3], N=update.sp.table[3,2])  #QD
  update.sp.table[4,4] <- get_stand(BA=update.sp.table[4,3], N=update.sp.table[4,2])  #QD
  update.sp.table[5,4] <- get_stand(BA=update.sp.table[5,3], N=update.sp.table[5,2])  #QD
  update.sp.table[1:5,5] <- NA

  return(list(Dd=Dd, Rec.SP=Rec.SP, sp.table=update.sp.table))

}

# Note: adjustment of VTHA is not included here.
