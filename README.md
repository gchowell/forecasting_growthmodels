# forecasting_growthmodels
 A Matlab toolbox for fitting and forecasting epidemic trajectories using phenomenological growth models


<p> It carries out the following tasks: </p> 
<ul>
    <li>fitting models to time series data,</li>
    <li>estimation of parameters with quantified uncertainty,</li>
    <li>plotting the best fit of the model and parameter estimates,</li>
    <li>plotting forecasts from the best fit model,</li>
    <li>generates and plots performance metrics of the forecasts,</li>
    <li>generates quantiles associated with the model fit and the forecast</li>

</ul>

<p> Additional features include:</p>

<ul>
    <li>fitting models using different parameter estimation approaches (least-squares, maximum likelihood estimation),</li>
    <li>fitting models using assuming different error structures (normal, Poisson, negagive binomial),</li>
    <li>user can select the underlying function for the sub-epidemic building block (generalized-growth model (GGM), generalized-logistic model (GLM), Richards model, Gompertz model),</li>
    <li> User can conduct multiple fits of the model to the data through a rolling-window analysis </li>
    
</ul>
    
# Installation requirements

The forecasting_growthmodels toolbox requires a MATLAB installation.

# Fitting the model to your data

To use the toolbox to fit a model to your data, you just need to:

<ul>
    <li>download the code </li>
    <li>create input folder where your time series data is located </li>
    <li>create output folder where the output files will be stored</li>   
    <li>open a MATLAB session </li>
    <li>define the model parameter values and time series parameters by editing <code>options_fit.m</code> </li>
    <li>run the function <code>Run_Fit_GrowthModels.m</code> </li>
</ul>
  
# Plotting the best model fit and parameter estimates

After fitting the model to your data, you can use the toolbox to plot the best model fit and parameter estimates as follows:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options_fit.m</code></li>
    <li>run the function <code>plotFit_GrowthModels.m</code> </li>
</ul>
    
# Plotting forecasts of the best fit model

After fitting the model to your data, you can use the toolbox to plot forecasts derived from the best fit model as follows:

<ul>
    <li>define the model parameter values and time series parameters by editing <code>options_forecast.m</code></li>
    <li>run the function <code>plotForecast_GrowthModels.m</code></li>
</ul>
    
# Publications

<ul>
    
<li> Chowell, G. (2017). Fitting dynamic models to epidemic outbreaks with quantified uncertainty: A primer for parameter uncertainty, identifiability, and forecasts. Infectious Disease Modelling, 2(3), 379-398. </li>

<li> Roosa, K., Lee, Y., Luo, R., Kirpich, A., Rothenberg, R., Hyman, J. M., ... & Chowell, G. (2020). Real-time forecasts of the COVID-19 epidemic in China from February 5th to February 24th, 2020. Infectious Disease Modelling, 5, 256-263.</li>

<li> Pell, B., Kuang, Y., Viboud, C., & Chowell, G. (2018). Using phenomenological models for forecasting the 2015 Ebola challenge. Epidemics, 22, 62-70. <li>

<li> Bürger, R., Chowell, G., & Lara-Díıaz, L. Y. (2019). Comparative analysis of phenomenological growth models applied to epidemic outbreaks. Mathematical Biosciences and Engineering, 16(5), 4250-4273. </li>

</ul>

# Disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.  
