# NOTHOPACK
A Growth and Yield Simulator for Nothofagus Second Growth Forests. Version 1.0.

## Instalation

In order to install **NOTHOPACK** you need:

````md
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("sgezan/Nothopack")
````

Once installed, to load this library in R you run the command:

````md
library(Nothopack)
````

***

### Description

**NOTHOPACK** is a growth model (G&Y) simulator for second-growth forests of *Nothofagus obliqua* (roble), *N. alpina* (raulí), and *N. dombeyi* (coihue), which are among the most important native mixed forests in Chile. They form the RORACO forest type that is present approximately between the 36 and 42 degrees S latitudes in both the Chilean Andes and the coastal mountain range with some fragments in Argentina. At the present, the RORACO forest type covers 1.96 million hectares, around 10% of the native forested area of Chile.

This simulator is the product of a large FONDEF project (D97I1065) that started in 1999 at the Universidad Austral under the direction of Dr. Alicia Ortega, financed by the Chilean Government. Under this 3-year project, several efforts were done to collect data for RORACO from previously studies. However, the largest effort was the development of two stratified sampling networks with temporary and permanent plots established across the complete range of this forest type in Chile. Efforts to build this simulator have continued over the years mainly at the Universidad Austral (Chile) and at the University of Florida (United States). Hence, **NOTHOPACK** is the result of combining all these resources and efforts into a single product, that has been put togehter into an R library to be used by foresters, and researchers to facilitate and improve the management and sustainability of this important Chilean forest type. 

This G&Y model has several properties, but the main ones are:

+ Inventory processing
+ Whole-stand level simulation
+ Individual-tree level simulation
+ Compatibility (stand-tree) simulation
+ Thinning simulation

## Authors

Salvador A. Gezan   <forestats.sg@gmail.com> \
Paulo C. Moreno     <paulo.moreno@ciep.cl>  \
Sebastian Palmas    <palmasforest@gmail.com>  \
Alicia Ortega       <aortega@uach.cl>

This R library and associated documentation can be accessed at: \
https://github.com/sgezan/Nothopack/

## Reference

**Gezan, S.A., Moreno, P.C., Palmas, S., and A. Ortega.** 2020. 
NOTHOPACK: A Growth and Yield Simulator for Nothofagus Second Growth Forests. Version 1.0.
Universidad Austral de Chile, Valdivia, Chile.
