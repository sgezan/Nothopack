#' Estimates the missing total height of trees from a given plot
#'
#' Estimates the total height of trees (m) that have missing height
#' that come from fitting a regression model on inventory data. For the missing
#' trees there are two methods: 1) to estimate tree heights according to a parametrized
#' height-DBH model, or 2) to estimate tree heights by fitting a simple height-DBH model
#' that requires at least 10 measurements. Missing values are indentified as 'NA', and
#' regression is fitted without making distinction of species.
#' All input vectors should correspond to the same tree ID.
#'
#' @param DBH Vector of diameter at breast height (cm)
#' @param HT Vector of total height (m) (with measured and/or missing values)
#' @param DOM.SP Dominant specie (1:Rauli, 2:Roble, 3:Coigue)
#' @param ZONE Growth zone (1, 2, 3, 4)
#' @param HD Dominant height (m) of dominant specie in the current stand
#' @param QD Quadratic diameter (cm) of the stand
#' @param method Method to use to estimate missing heights. 1:Parametrized DBH-height model
#'          (which requires input: DBH, DOM.SP, ZONE, HD, QD), and 2:Fits a simple DBH-height model from
#'          available measurements using the equation: \eqn{ln(HT) = b0 + b1/DBH} (default = 1)
#'
#' @return A list containing the following:
#' \itemize{
#' \item \code{HTFIN} A vector of final tree heigths (m), replacing the missing values for estimated heigths and
#'                    retaining the observed heights
#' \item \code{r2} A value with the coefficient of determination from the fitting the DBH-height model when method = 2
#' }
#'
#' @author
#' P. Someda-Dias, and S.A. Gezan
#'
#' @examples
#' # Example 1: Method 1 - Parametrized DBH-height model
#' DBH <- c(9.3, 11.1, 15.5, 9, 14.8, 27.3, 11.4, 6.6, 12.6, 17.5,
#'          6.3, 7.2, 11.5, 13.6, 7.3, 12, 11.9, 8.1, 7.6, 5)
#' HT <- c(11.8, 12.3, NA, NA, 15.3, 18, 12, NA, 14.5, NA, NA,
#'         NA, NA, NA, 10.3, 14.6, NA, NA, NA, NA)
#' (ht.est <- tree.ht(DBH=DBH, HT=HT, DOM.SP=2, ZONE=2, HD=15, QD=12, method=1)$HTFIN)
#' plot(DBH, ht.est, ylab='HT (m)', xlab='DBH (cm)', xlim=c(0,30), ylim=c(0,20),
#'      type='p', col='blue')
#'
#' # Example 2: Method 2 - Simple DBH-height model
#' DBH <- c(9.3, 11.1, 15.5, 9, 14.8, 27.3, 11.4, 6.6, 12.6, 17.5,
#'          6.3, 7.2, 11.5, 13.6, 7.3, 12, 11.9, 8.1, 7.6, 5)
#' HT <- c(11.8, 12.3, NA, NA, 15.3, 18, 12, NA, 14.5, NA, NA,
#'         NA, NA, NA, 10.3, 14.6, NA, NA, NA, NA)
#' (ht.est <- tree.ht(DBH=DBH, HT=HT, DOM.SP=2, ZONE=2, HD=15, QD=12, method=2)$HTFIN)
#' (r2 <- tree.ht(DBH=DBH, HT=HT, DOM.SP=2, ZONE=2, HD=15, QD=12, method=2)$r2)
#' plot(DBH, ht.est, ylab='HT (m)', xlab='DBH (cm)',
#'      xlim=c(0,30), ylim=c(0,20), type='p', col='darkblue')

tree.ht  <-  function(DBH, HT, DOM.SP=NA, ZONE=NA, HD=NA, QD=NA, method=1){

  if(length(DBH)==length(HT) & sum(is.na(DBH))==0){

    if(method==1){

      if(is.na(DOM.SP)==T || is.na(ZONE)==T || is.na(HD)==T || is.na(QD)==T){
        stop("Warning - Incomplete information for tree_ht. Please check input.")
      }

      aux <- data.frame(DBH=DBH, HT=HT, HTEST=NA, HTFIN=NA)
      aux$HTEST <- height_param(DOM.SP=DOM.SP, ZONE=ZONE, HD=HD, QD=QD, DBH=DBH)
      aux$HTFIN <- ifelse(is.na(aux$HT),aux$HTEST,aux$HT)
      r2 <- NA

    }

    if(method==2){

      aux <- data.frame(DBH=DBH, HT=HT)
      aux2 <-subset(aux, HT!='NA')
      model <- stats::lm(log(aux2$HT)~I(1/aux2$DBH),na.action='na.omit')
      r2 <- summary(model)$r.squared
      HTEST <- exp(model$coefficients[1]+model$coefficients[2]/DBH)
      aux$HTFIN <- ifelse(is.na(aux$HT),HTEST,aux$HT)

    }

  } else {

    stop("Warning - Incomplete information for tree_ht. Please check input.")

  }

  return(list(HTFIN=aux$HTFIN, r2=r2))

}
