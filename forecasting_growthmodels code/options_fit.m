% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
function [cadfilename1,caddisease,datatype, dist1, numstartpoints,M,flag1,model_name1,fixI0,printscreen1,windowsize1,tstart1,tend1]=options_fit

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% The time series data file contains the incidence curve of the epidemic of interest. 
% The first column corresponds to time index: 0,1,2, ... and the second
% column corresponds to the temporal incidence data. If the time series file contains cumulative incidence count data, 
% the name of the time series data file must start with "cumulative".

cadfilename1='Most_Recent_Timeseries_US-CDC'; % Name of the data file containing the incidence curve

caddisease='monkeypox'; % string indicating the name of the disease related to the time series data

datatype='cases'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=0; % Type of estimation method

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % Define dist1 which is the type of error structure. See below:

%dist1=0; % Normal distribution to model error structure (method1=0)
%dist1=1; % Poisson error structure (method1=0 OR method1=1)
%dist1=2; % Neg. binomial error structure where var = factor1*mean where
                  % factor1 is empirically estimated from the time series
                  % data (method1=0)
%dist1=3; % MLE (Neg Binomial) with VAR=mean+alpha*mean  (method1=3)
%dist1=4; % MLE (Neg Binomial) with VAR=mean+alpha*mean^2 (method1=4)
%dist1=5; % MLE (Neg Binomial)with VAR=mean+alpha*mean^d (method1=5)


switch method1
    case 1
        dist1=1;
    case 3
        dist1=3;
    case 4
        dist1=4;
    case 5
        dist1=5;
end


numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

M=200; % number of bootstrap realizations to characterize parameter uncertainty

% <==============================================================================>
% <============================== Growth model =====================================>
% <==============================================================================>

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=GLM; % Growth model considered in the epidemic trajectory

model_name1='GLM';  % A string provided by the user for the name of the model

% <==================================================================================>
% <=============================== other parameters=======================================>
% <==================================================================================>

fixI0=1; % 0=Estimate the initial number of cases; 1 = Fix the initial number of cases according to the first data point in the time series

printscreen1=1;  % flag (1/0) to indicate whether we want to print plots with the results

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=8;  % moving window size

tstart1=1; % time point for the start of rolling window analysis

tend1=1;  %time point for the end of the rolling window analysis
