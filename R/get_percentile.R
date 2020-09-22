#' Calculates percentiles for DBH based on vectors of DBH and N
#'
#' If calculates percentiles for DBH based on the frequency
#' table provided by vectors of DBH and number of trees (N). Used for
#' expansion factors (or N) corrections.
#'
#' @param DBH Vector of diameter breast height values
#' @param FT Vector of expansion factor values (or N, in case of diameter distribution)
#' @param percentiles list of desired percentile values
#'
#' @return Data frame with desired DBH percentile values
#'
#' @author
#' S. Palmas, S.A. Gezan and P. Moreno
#'
#' @examples
#' DBH <- c(36.5, 2.8, 0, 2.4)
#' N <- c(464, 23, 0, 48)
#' get_percentile(DBH=DBH, FT=N, percentiles=c(0.15, 0.85))

get_percentile <- function(DBH=NA, FT=NA, percentiles=c(0.15, 0.85)){

    df <- data.frame(DBH=DBH, FT=FT)
    df <- df[df$FT>0,]        # Eliminate any line that is 0 (anywhere)
    df <- rbind(df,c(0.0))    # Adding the initial row of 0
    df <- df[order(df$DBH),]  # Sorting from small diameter to largest
    df$CumFT <- cumsum(df$FT) # Adding cummulative sum

    Cperc <- sum(df$FT, na.rm = TRUE) * percentiles    #Percentiles adjusted for total FT
    places <- as.numeric(cut(Cperc, breaks = df$CumFT, right = FALSE)) #Places of the percentiles inside the CumFT array
    values <- (df$CumFT[places]*df$DBH[places]+df$DBH[places+1]*(Cperc - df$CumFT[places]))/Cperc    #Calculating the adjusted percentiles

    dfin <- data.frame(Percentile=percentiles, DBH=values)

  return(dfin)
}
