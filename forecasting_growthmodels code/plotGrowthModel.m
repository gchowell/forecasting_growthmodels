function plotGrowthModel(flag1_pass,params0_pass,I0_pass,windowsize1_pass)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% Fitting and forecasting model to epidemic data with quantified uncertainty

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method

% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

[cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP,flag1_INP,model_name1_INP,fixI0_INP,windowsize1_INP,tstart1_INP,tend1_INP]=options_fit;

% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

DT=1;

caddisease=caddisease_INP;

datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

%method1=0; % Type of estimation method: 0 = LSQ

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squard=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;


numstartpoints=numstartpoints_INP; % Number of initial guesses for optimization procedure using MultiStart

M=M_INP; % number of bootstrap realizations to characterize parameter uncertainty


% <==============================================================================>
% <============================== Growth model =====================================>
% <==============================================================================>

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

model_name1=model_name1_INP;

% <==================================================================================>
% <=============================== other parameters=======================================>
% <==================================================================================>

fixI0=fixI0_INP; % 0=Estimate the initial number of cases; 1 = Fix the initial number of cases according to the first data point in the time series

% <==============================================================================>
% <======================== Load epidemic data ========================================>
% <==============================================================================>

% data=load(strcat('./input/',cadfilename1,'.txt'));
% 
% if strcmp('CUMULATIVE',upper(cadfilename1(1:10)))==1
% 
%     data(:,2)=[data(1,2);diff(data(:,2))]; % Incidence curve
% 
% end

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=1; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=0; %forecast horizon (number of data points ahead)

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

if exist('flag1_pass','var')==1 & isempty(flag1_pass)==0
    flag1=flag1_pass;
else
    flag1=flag1_INP;
end

params0=params0_pass;

I0=I0_pass;

if exist('windowsize1_pass','var')==1 & isempty(windowsize1_pass)==0
    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end

%params0=initialParams(data(:,2),flag1);

r=params0(1);
p=params0(2);
a=params0(3);
K=params0(4);

timevect=0:1:windowsize1;

options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);

[~,F]=ode15s(@modifiedLogisticGrowth,timevect,I0,options,r,p,a,K,flag1);

fitcurve=abs([F(1,1);diff(F(:,1))]);

plot(timevect,fitcurve,'b-');





