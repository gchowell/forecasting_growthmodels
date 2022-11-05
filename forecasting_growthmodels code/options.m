
function [cadfilename1,DT, dist1, numstartpoints,M,flag1,model_name1,fixI0,getperformance,forecastingperiod, printscreen1,windowsize1,tstart1,tend1]=options_fit

% Parameters for FittingModelToData.m function that fits growth models to
% time series data with quantified uncertainty


% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1='sim-GLM-poiss';  %name of the time series data file

DT=1;  %temporal resolution of the data (DT=1 for daily data).

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=0; % Type of estimation method: 0 = LSQ

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squared=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % Define dist1 which is the type of error structure:

%dist1=0; % Normnal distribution to model error structure
%dist1=1; % error structure type (Poisson=1; NB=2)
%dist1=3; % VAR=mean+alpha*mean;
%dist1=4; % VAR=mean+alpha*mean^2;
%dist1=5; % VAR=mean+alpha*mean^d;

numstartpoints=6; % Number of initial guesses for optimization procedure using MultiStart

M=100; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== Growth model =====================================>
% <==============================================================================>

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=GLM; % Sequence of subepidemic growth models considered in epidemic trajectory

model_name1='GLM';  % name of the model

% <==================================================================================>
% <=============================== other parameters=======================================>
% <==================================================================================>

fixI0=1; % 0=Estimate the initial number of cases; 1 = Fix the initial number of cases according to the first data point in the time series


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=30; %forecast horizon (number of data points ahead)

printscreen1=1;  % print plots with the results

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=20;  %moving window size
tstart1=10; % time of start of rolling window analysis
tend1=tstart1+windowsize1-1;  %time end of the rolling window analysis
%tend1=length(data(:,1));

