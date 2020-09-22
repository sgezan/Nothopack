#' Calculates tree volume (inside bark) based on taper equations for different products according to specifications
#'
#' Individual tree volume for different products, according to a table of product specifications,
#' is calculated over all individuals provided in a tree list (based on DBH, HT, SPECIES).
#' A new tree list with additional columns is generated, which can be used to obtain a stand table.
#' Note: the table of product specificiacionts needs to be ordered from most to least important.
#'
#' @param tree.list Data frame containing the following vectors:
#'          ID (tree id), SPECIES (tree specie), DBH (diameter breast height, cm),
#'          HT (total tree height, m). Additional columns can be included but these will be ignored
#' @param ZONE Growth zone of the corresponding stand
#' @param prod.list Data frame containing the following vectors:
#'          Product (product abbreviation), Desc (description of product),
#'          Dmin (minimum diameter for product, cm), Length (length of section, m)
#' @param T_model Model selected for taper equation (1:zone specific, 2:all zones) (default = 2)
#' @param stump Length of stump to discount (default = 0.3 m)
#'
#' @return Updated tree list with additional columns of individual tree
#' volume products (m3) (one column per product)
#'
#' @author
#' S.A. Gezan
#'
#' @references
#' Gezan, S.A. and Moreno, P. (2000b). Informe Modelos de Volumen.
#' Reporte Interno. Projecto FONDEF D97I1065. Chile
#'
#' @examples
#' # Example: Calculation of volumes for 5 different products with their stand table
#' tree.data <- plot_example
#' head(tree.data,10)  # Original data
#' tree.inv <- input_module(ZONE=1, AD=28, type='tree', area=500, tree.list=tree.data, T_model=2)
#' head(tree.inv$tree.list, 10)
#' Product <- c('SL', 'SS', 'PP', 'IF', 'RES')
#' Desc <- c('Saw-log (large)','Saw-log (short)','Panel','Industrial Firewood','Residual')
#' Dmin <- c(18, 16, 10, 5, 0)
#' Length <- c(4.1, 3.2, 2.44, 1.22, 99)
#' (Vprod <- data.frame(Product,Desc,Dmin,Length))
#' updated.tree.list <- VPRODmodule(tree.list=tree.inv$tree.list, ZONE=tree.inv$ZONE,
#'                                prod.list=Vprod, T_model=2)
#' head(updated.tree.list, 10)
#' tree.inv$sp.table
#' plot.extra <- stand_parameters(plotdata=updated.tree.list, area=500,
#'                   extra=c(8,9,10,11,12,13,14))
#' plot.extra$sd

VPRODmodule <- function(tree.list=NA, ZONE=NA, prod.list=NA, T_model=2, stump=0.3){

  if(all(is.na(prod.list))){
    stop('A product list needs to be provided')
  }

  products <- prod.list
  n.p<-length(products$Product)
  n.l<-length(tree.list$ID)

  # Volume first product (later more generic)
  VOL.p <- matrix(data=0,nrow=n.l,ncol=n.p+1)
  #head(VOL.p)

  for (i in 1:n.l) {

    SPECIESi<-tree.list$SPECIES[i]
    DBHi<-tree.list$DBH[i]
    HTi<-tree.list$HT[i]

    stump0 <- stump
    p <- 1
    pass <- 0
    while (pass==0) {
      hi <- stump0 + products$Length[p]
      di <- get_taper(SPECIES=SPECIESi, ZONE=ZONE, DBH=DBHi, HT=HTi, model=T_model, hi=hi)$di
      if (di >= products$Dmin[p] & hi <= HTi) {
        #pass <- 1
        # Example 1: Calculates tree volume for diameter limit of 5 cm (with stump of 0.3 m)
        VOL.p[i,p] <- VOL.p[i,p] + Vmodule_individual(SPECIES=SPECIESi, ZONE=ZONE,
                                                      DBH=DBHi, HT=HTi, dmin=di, model=T_model, stump=stump0)
        stump0 <- hi
        #print(stump0)
      } else {
        p <- p+1
        if (p > n.p) {
          pass<-1
          VOL.p[i,p-1] <- VOL.p[i,p-1] + Vmodule_individual(SPECIES=SPECIESi, ZONE=ZONE,
                                              DBH=DBHi, HT=HTi, dmin=0, model=T_model, stump=stump0)
        }
      }
    }
    VOL.p[i,p] <-sum(VOL.p[i,])
  }

  nam <- paste0(rep('VOLP_',p-1), products$Product, sep='')
  nam[p] <- "VTOT"
  colnames(VOL.p) <- nam
  #head(VOL.p)

  updated.list <- cbind(tree.list,VOL.p)

  return(updated.list)
}
