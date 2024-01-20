# forecasting_growthmodels
 A Matlab toolbox for fitting and forecasting epidemic trajectories using phenomenological growth models


<p> It carries out the following tasks: </p> 
<ul>
    <li>fitting models to time series data,</li>
    <li>estimation of model parameters, Monte Carlo standard errors (MCSES) of the parameter estimates, doubling times, and reproduction numbers with quantified uncertainty,</li>
    <li>plotting the best fit of the model and parameter estimates,</li>
    <li>plotting forecasts from the best-fit model,</li>
    <li>generates and plots performance metrics of the forecasts,</li>
    <li>generates quantiles associated with the model fit and the forecast</li>
    <li>conducts rolling window analyses of parameter estimates during specific time periods and window sizes</li>
    
</ul>

<p> Additional features include:</p>

<ul>
    <li>fitting models using different parameter estimation approaches (least-squares, maximum likelihood estimation),</li>
    <li>fitting models using assuming different error structures (normal, Poisson, negative binomial),</li>
    <li>user can select among different growth models including the generalized-growth model (GGM), generalized-logistic growth model (GLM), Gompertz model, Richards model, and the generalized Richards model (GRM),</li>
    <li> User can conduct multiple fits of the model to the data through a rolling-window analysis </li>
    
</ul>
    
# Installation requirements

The forecasting_growthmodels toolbox requires a MATLAB installation.

# Fitting the model to your data

To use the toolbox to fit a model to your data, you just need to:

<ul>
    <li>download the code </li>
    <li>create 'input' folder in your working directory where your data is located </li>
    <li>create 'output' folder in your working directory where the output files will be stored</li>   
    <li>open a MATLAB session </li>
    <li>define the model parameter values and time series parameters by editing <code>options_fit.m</code> </li>
    <li>run the function <code>Run_Fit_GrowthModels.m</code> </li>
</ul>
  
# Plotting the best model fit and parameter estimates

After fitting the model to your data using the function <code>Run_Fit_GrowthModels.m</code>, you can use the toolbox to plot the best model fit and parameter estimates as follows:

<ul>
    <li>run the function <code>plotFit_GrowthModels.m</code> </li>
</ul>

The function also outputs files with parameter estimates, the best fit of the model, and the performance metrics for the calibration period.

# Generating forecasts

To use the toolbox to fit a model to your data and generate a forecast, you just need to:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options_forecast.m</code> </li>
    <li>run the function <code>Run_Forecasting_GrowthModels.m</code> </li>
</ul>
  
# Plotting forecasts and performance metrics based on the best-fit model

After running <code>Run_Forecasting_GrowthModels.m</code>, you can use the toolbox to plot forecasts and performance metrics derived from the best fit model as follows:

<ul>
    <li>run the function <code>plotForecast_GrowthModels.m</code></li>
</ul>

The function also outputs files with parameter estimates, the fit and forecast of the model, and the performance metrics for the calibration period and forecasting periods.

# Publications

<ul>

<li>Chowell, G., Bleichrodt, A., Dahal, S. et al. GrowthPredict: A toolbox and tutorial-based primer for fitting and forecasting growth trajectories using phenomenological growth models. Sci Rep 14, 1630 (2024). https://doi.org/10.1038/s41598-024-51852-8 </li>
    
<li> Chowell, G. (2017). Fitting dynamic models to epidemic outbreaks with quantified uncertainty: A primer for parameter uncertainty, identifiability, and forecasts. Infectious Disease Modelling, 2(3), 379-398. </li>

<li> Bürger, R., Chowell, G., & Lara-Díıaz, L. Y. (2019). Comparative analysis of phenomenological growth models applied to epidemic outbreaks. Mathematical Biosciences and Engineering, 16(5), 4250-4273. </li>

</ul>

# Disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.  
