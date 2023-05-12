% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
function [cadfilename1,caddisease,datatype, dist1, numstartpoints,B,flag1,model_name1,fixI0,windowsize1,tstart1,tend1]=options_fit

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
% Located in the input folder, the time series data file is a text file with extension *.txt. 
% The time series data file contains the incidence curve of the epidemic of interest. 
% The first column corresponds to time index: 0,1,2, ... and the second
% column corresponds to the temporal incidence data. If the time series file contains cumulative incidence count data, 
% the name of the time series data file must start with "cumulative".

cadfilename1='Most_Recent_Timeseries_US-CDC'; % String variable indicating the name of the data file containing the time-series data.

caddisease='monkeypox'; % string variable indicating the name of the disease or subject related to the time series data

datatype='cases'; % string variable indicating the nature of the data (cases, deaths, hospitalizations, etc)

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=0; % This integer variable indicates the parameter estimation method employed to estimate the parameters from data. 
% The following estimation methods are available:

% Nonlinear least squares (LSQ)=0,
% MLE Poisson=1,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

dist1=0; % This integer variable indicates the error structure assumptions. The following error structure assumptions are available:

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

numstartpoints=10; % This variable defines the number of different initial guesses for the optimization procedure using Multistart 
% in its search for the globally optimal set of parameters.

B=300; % Number of bootstrap realizations utilized to characterize parameter uncertainty.

% <==============================================================================>
% <============================== Growth model =====================================>
% <==============================================================================>

EXP=-1;  % -1 = EXP
GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards
GOM=5; % 5 = Gompertz


flag1=GLM; % Integer variable indicating the growth model that will be fit to the time-series data.

model_name1='GLM';  % A string variable indicating the name of the model.

fixI0=0; % Boolean variable indicating whether initial value in the time-series will be estimated or fix according to the first data point in the time series.

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=31;  % Integer variable indicating the moving window size

tstart1=1; % Integer variable indicating the time point for the start of rolling window analysis

tend1=2;  %Integer variable indicating the time point for the end of the rolling window analysis
