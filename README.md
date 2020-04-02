# COVID19
An attempt to model the COVID19 data from the [Johns Hopkins Github Repository](https://github.com/CSSEGISandData) by a SIRD model for various countries such as Germany, Switzerland and the US. This Github repository contains a shiny app that estimates the model parameters and then simulates how to disease would progress based on these estimates.

A live version of the shiny app can be found [here](https://joelmeili.shinyapps.io/COVID19).

## Features of the shiny app
- estimation of the SIRD model coefficients with adjustable estimation interval and forecasting horizon
- multiple plots that illustrate the spread of the disease for various countries
- multiple plots that illustrate the estimated SIRD curves for various countries

## SIRD model
![](sird.jpg)

## SIRD model equations
![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdS%7D%7Bdt%7D%20%3D%20-%5Cbeta%20I%28t%29%20S%28t%29) <br/>
![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdI%7D%7Bdt%7D%20%3D%20%5Cbeta%20I%28t%29%20S%28t%29%20-%20%5Cgamma%20I%28t%29%20-%20%5Cdelta%20I%28t%29) <br/>
![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdR%7D%7Bdt%7D%20%3D%20%5Cgamma%20I%28t%29) <br/>
![](http://latex.codecogs.com/gif.latex?%5Cfrac%7BdD%7D%7Bdt%7D%20%3D%20%5Cdelta%20I%28t%29) <br/>
![](http://latex.codecogs.com/gif.latex?R_0%20%3D%20N%20%5Cfrac%7B%5Cbeta%7D%7B%5Cgamma%20&plus;%20%5Cdelta%7D)

where N is the size of the population which is assumed to be constant.

## Model estimation
To estimate the model parameters a RMSE estimator was used that compares the simulated number and the true number for infected, recovered and dead cases.