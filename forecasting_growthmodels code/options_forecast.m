
% options_forecast.m — GrowthPredict
% Put this file in your working directory (same level as ./input and ./output).
% The GrowthPredict user function Run_Forecasting_GrowthModels.m reads these
% variables from the workspace.

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, flag1, model_name1, fixI0, getperformance, forecastingperiod, windowsize1, tstart1, tend1] = options_forecast

% OPTIONS_FORECAST  GrowthPredict options for forecasting with growth models
%
% Overview
%   Returns configuration to (i) calibrate the chosen growth model on rolling
%   windows and (ii) generate out-of-sample forecasts with quantified uncertainty.
%
% Usage
%   [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, ...
%    flag1, model_name1, fixI0, getperformance, forecastingperiod, ...
%    windowsize1, tstart1, tend1] = options_forecast;
%
% Returns
%   cadfilename1      char     Base name of the input file in ./input ('<cadfilename1>.txt')
%   caddisease        char     Disease/subject label for outputs
%   datatype          char     Data tag (e.g., 'cases' | 'deaths' | 'hospitalizations')
%   dist1             int      Error model (see mapping below)
%   numstartpoints    int      MultiStart initial points for global search
%   B                 int      Bootstrap replicates for uncertainty
%   flag1             int      Growth-model code (see “Growth model choices”)
%   model_name1       char     Human-readable model name matching flag1
%   fixI0             double   1 = fix initial observed value to first datum; 0 = estimate it
%   getperformance    logical  1=compute forecast performance metrics; 0=skip
%   forecastingperiod int      Forecast horizon (steps ahead)
%   windowsize1       int      Rolling-window length (time steps)
%   tstart1           int      Start index of the first rolling window
%   tstart1         int        Start index of the first rolling window

%
% Input data (./input)
%   Text file '<cadfilename1>.txt' with two columns, NO header:
%     Col 1: time index  (0,1,2,...)
%     Col 2: observed incidence (nonnegative)

%   If the series is cumulative, the filename must BEGIN with the word
%   'cumulative' (case-insensitive; e.g., 'cumulative-mpox.txt' or
%   'CumulativeCases.txt'). Such files are converted to incidence via
%   [C(1); diff(C)] on load. Any other name is treated as incidence as-is,
%   so a cumulative series without this prefix will be silently misread.
%%
% Estimation & error models
%   The global 'method1' selects the estimator; 'dist1' sets/weights the observation model.
%     method1 = 0  LSQ with dist1 ∈ {0,1,2} as weighting (Normal / Poisson-like / NegBin-like)
%     method1 = 1  MLE Poisson                          → dist1 := 1 (automatic)
%     method1 = 3  MLE NegBin: var = mean + α·mean      → dist1 := 3 (automatic)
%     method1 = 4  MLE NegBin: var = mean + α·mean^2    → dist1 := 4 (automatic)
%     method1 = 5  MLE NegBin: var = mean + α·mean^d    → dist1 := 5 (automatic)
%
% Growth model choices (flag1)
%   EXP=-1 (exponential), GGM=0 (generalized growth),
%   GLM=1 (generalized logistic), GRM=2 (generalized Richards),
%   LM=3 (logistic), RICH=4 (Richards), GOM=5 (Gompertz).
%   Set model_name1 to the matching abbreviation (e.g., 'GLM').
%
% Notes
%   • Choose forecastingperiod to match your application (e.g., 4 weeks if weekly data).
%   • getperformance=1 writes forecast accuracy metrics to ./output (if implemented).
%   • Keep file/disease labels ASCII if you need cross-platform filename compatibility.


% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Global variable used to define the parameter estimation method.
global method1; % Parameter estimation method

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% The time series data file contains the incidence curve of interest (e.g., new cases per unit of time).
% - The first column corresponds to the time index (e.g., 0, 1, 2, ...).
% - The second column contains the temporal incidence data.

% Note: If the time series file contains cumulative count data, its name must
%       begin with the word "cumulative" (case-insensitive; a hyphen after it
%       is conventional but not required). Files matching this prefix are
%       differenced into incidence on load; all other files are used as-is.

cadfilename1 = 'Most_Recent_Timeseries_US-CDC'; % Name of the time-series data file
caddisease = 'Mpox';                            % Name of the disease or subject related to the data
datatype = 'cases';                             % Type of data (e.g., cases, deaths, hospitalizations)

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Estimation method options:
% 0 - Nonlinear least squares (LSQ)
% 1 - Maximum Likelihood Estimation (MLE) Poisson
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

method1 = 0; % Default estimation method: Nonlinear least squares (LSQ)

% Error structure options based on method1:
% 0 - Normal distribution (method1 = 0)
% 1 - Poisson error structure (method1 = 0 or 1)
% 2 - Negative Binomial (VAR = factor1 * mean, empirically estimated)
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

dist1 = 0; % Default error structure: Normal distribution
switch method1
    case 1
        dist1 = 1; % Poisson error structure
    case 3
        dist1 = 3; % Negative Binomial (VAR = mean + alpha*mean)
    case 4
        dist1 = 4; % Negative Binomial (VAR = mean + alpha*mean^2)
    case 5
        dist1 = 5; % Negative Binomial (VAR = mean + alpha*mean^d)
end

% Optimization settings:
numstartpoints = 10; % Number of initial guesses for global optimization (Multistart)
B = 100;             % Number of bootstrap realizations for parameter uncertainty characterization

% <============================================================================>
% <========================== Growth Model ===================================>
% <============================================================================>
% Growth-model definitions for the cumulative trajectory C(t):
% -1: Exponential Growth Model (EXP)
%     dC/dt = r*C
%  0: Generalized Growth Model (GGM)
%     dC/dt = r*C^p
%  1: Generalized Logistic Growth Model (GLM)
%     dC/dt = r*C^p*(1 - C/K)
%  2: Generalized Richards Model (GRM)
%     dC/dt = r*C^p*(1 - (C/K)^a)
%  3: Logistic Growth Model (LM)
%     dC/dt = r*C*(1 - C/K)
%  4: Richards Model (RICH)
%     dC/dt = r*C*(1 - (C/K)^a)
%  5: Gompertz Model (GOM), time-dependent growth-rate form
%     dC/dt = r*C*exp(-a*t); K is not used in the current implementation
%
% The fitted/forecast incidence vector is constructed as [C(t_1); diff(C(t))].

EXP = -1;  GGM = 0;  GLM = 1;  GRM = 2;  LM = 3;  RICH = 4;  GOM = 5;

flag1 = GLM;         % Selected growth model: Generalized Logistic Growth Model (GLM)
model_name1 = 'GLM'; % Name of the selected model
fixI0 = 1;           % Boolean: Fix the initial value to the first data point (true) or estimate it (false)

% <============================================================================>
% <====================== Forecasting Parameters =============================>
% <============================================================================>
% Parameters for forecasting analysis:
% - getperformance: Boolean to enable/disable forecasting performance metrics
% - forecastingperiod: Time horizon for forecasting (number of time units ahead)

getperformance = 1;    % Enable forecasting performance metrics (1 = yes, 0 = no)
forecastingperiod = 4; % Forecast horizon: Number of time units ahead

% <============================================================================>
% <======= Parameters for Rolling Window Analysis ===========================>
% <============================================================================>
% Parameters for rolling window analysis:
% - windowsize1: Size of the moving window
% - tstart1: Start time point for rolling window analysis
% - tend1: End time point for rolling window analysis

windowsize1 = 10; % Size of the rolling window (e.g., 10 time units)
tstart1 = 1;      % Start time point for rolling window analysis
tend1 = 1;        % End time point for rolling window analysis

end
