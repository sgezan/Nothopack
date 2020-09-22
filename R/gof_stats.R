#' Calculates goodness-of-fit statistics
#'
#' It obtains goodness-of-fit statistics (r2emp, RMSE, BIAS, RMSE%, BIAS%, Theil)
#' using provided observed and fitted values. Missing values are removed from the list.
#'
#' @param obs vector of observed values
#' @param pred vector of predicted values
#'
#' @return A table with goodness-of-fit statistics:
#' r2emp: empirical coefficient of correlation
#' RMSE (and RMSE%) root mean square error (raw and in percentage of mean of observed data)
#' BIAS (and BIAS%) bias (raw and in percentage of mean of observed data)
#' Theil's Inequality Coefficient (Ahlburg 1984)
#'
#' @references 
#' Alburhg, D.A. (1984). Forecast Evaluation and Improvement using Thiel's 
#' Decomposition. Journal of Forecasting 3(3):345-351
#' 
#' @author 
#' P. Moreno & S.A. Gezan
#' 
#' @examples
#' x <- c(1:20)
#' yobs <- 10 + 0.1*x + rnorm(20, mean=0, sd=0.2)
#' fit <- lm(yobs ~ x)
#' gof_stats(yobs, predict(fit))

gof_stats <- function(obs, pred){
  
  n<-length(obs)
  mean_obs  <- mean(obs, na.rm=TRUE)
  fit_r2emp <- round(1 - sum((obs - pred)^2, na.rm=TRUE) / sum((obs - mean_obs)^2, na.rm=TRUE), digits=3)
  fit_rmse  <- round(sqrt( sum( (obs - pred)^2, na.rm=TRUE) / (length(obs) - 1) ), digits=3)  # This was fixed
  fit_rmsep <- round(100 * fit_rmse/mean_obs, digits=3)
  fit_bias  <- round( sum( obs - pred, na.rm=TRUE)/ (length(obs)) , digits=3)
  fit_biasp <- round(100 * fit_bias/mean_obs, digits=3)
  fit_Theils<- round(sqrt((sum((obs - pred)^2, na.rm=TRUE)/sum((obs)^2, na.rm=TRUE))), digits=3)

  tabla <- matrix(c(n,fit_r2emp, fit_rmse, fit_rmsep, fit_bias, fit_biasp, fit_Theils),
                  ncol = 7)
  tabla <- as.data.frame(tabla)
  colnames(tabla) <- c('n','r2emp','RMSE','RMSE%','BIAS','BIAS%','Theil')

  return (tabla)
}

