
% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>
function [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, flag1, model_name1, fixI0, windowsize1, tstart1, tend1] = options_fit

% OPTIONS_FIT  GrowthPredict options for fitting single growth models
%
% Overview
%   Returns all configuration needed to calibrate a chosen growth model to a
%   univariate time series and quantify uncertainty (via parametric bootstrap).
%
% Usage
%   [cadfilename1, caddisease, datatype, dist1, numstartpoints, B, ...
%    flag1, model_name1, fixI0, windowsize1, tstart1, tend1] = options_fit;
%
% Returns
%   cadfilename1    char     Base name of the input file in ./input (expects '<cadfilename1>.txt')
%   caddisease      char     Disease/subject label for outputs (e.g., 'Mpox')
%   datatype        char     Data tag (e.g., 'cases' | 'deaths' | 'hospitalizations')
%   dist1           int      Error model (see “Estimation & error models” below)
%   numstartpoints  int      MultiStart initial points for global search
%   B               int      Bootstrap replicates for parameter uncertainty
%   flag1           int      Growth-model code (see “Growth model choices”)
%   model_name1     char     Human-readable model name matching flag1 (e.g., 'GLM')
%   fixI0           double   1 = fix initial observed value to first datum; 0 = estimate it
%   windowsize1     int      Rolling-window length (time steps)
%   tstart1         int      Start index of the first rolling window
%   tend1           int      Start index of the last rolling window
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
%
% Estimation & error models
%   The global 'method1' selects the estimator; 'dist1' sets/weights the observation model.
%     method1 = 0  Nonlinear least squares (LSQ)
%         choose dist1 ∈ {0,1,2}:
%           0 = Normal (homoscedastic LSQ)
%           1 = Poisson-like weighting (var ≈ mean; LSQ variant)
%           2 = NegBin-like weighting with var = factor1·mean (factor1 estimated empirically)
%     method1 = 1  MLE Poisson                          → dist1 := 1 (automatic)
%     method1 = 3  MLE NegBin: var = mean + α·mean      → dist1 := 3 (automatic)
%     method1 = 4  MLE NegBin: var = mean + α·mean^2    → dist1 := 4 (automatic)
%     method1 = 5  MLE NegBin: var = mean + α·mean^d    → dist1 := 5 (automatic)
%
% Growth model choices (flag1)
%   EXP = -1 (Exponential Growth), GGM = 0 (Generalized Growth),
%   GLM = 1 (Generalized Logistic Growth), GRM = 2 (Generalized Richards),
%   LM = 3 (Logistic Growth), RICH = 4 (Richards), GOM = 5 (Gompertz).
%
% Notes
%   • MultiStart (numstartpoints) helps avoid local minima for nonlinear models.
%   • Set fixI0=1 to anchor the initial observed state to the first data point.
%   • Rolling windows use indices in the time index, not calendar dates.


% <============================================================================>
% <=================== Declare Global Variables ==============================>
% <============================================================================>
% Global variables used throughout the function.
global method1; % Parameter estimation method

% <============================================================================>
% <========================= Dataset Properties ==============================>
% <============================================================================>
% The time series data file is a text file (*.txt) located in the input folder. 
% This file contains the incidence curve of interest (e.g., new cases per unit of time).
% - The first column corresponds to the time index: 0, 1, 2, ...
% - The second column contains the temporal incidence data.
% Note: If the time series file contains cumulative count data, its name must
%       begin with the word "cumulative" (case-insensitive; a hyphen after it
%       is conventional but not required). Files matching this prefix are
%       differenced into incidence on load; all other files are used as-is.

cadfilename1 = 'Most_Recent_Timeseries_US-CDC'; % Name of the data file containing the time-series data.
caddisease = 'Mpox';                            % Name of the disease or subject related to the time series.
datatype = 'cases';                             % Nature of the data (e.g., cases, deaths, hospitalizations).

% <============================================================================>
% <======================= Parameter Estimation ==============================>
% <============================================================================>
% Method used for parameter estimation:
% 0 - Nonlinear least squares (LSQ)
% 1 - Maximum Likelihood Estimation (MLE) Poisson
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

method1 = 0; % Default estimation method: Nonlinear least squares (LSQ).

% Error structure assumptions:
% 0 - Normal distribution (for method1 = 0)
% 1 - Poisson error structure (for method1 = 0 or 1)
% 2 - Negative Binomial (VAR = factor1 * mean, empirically estimated)
% 3 - MLE Negative Binomial (VAR = mean + alpha*mean)
% 4 - MLE Negative Binomial (VAR = mean + alpha*mean^2)
% 5 - MLE Negative Binomial (VAR = mean + alpha*mean^d)

dist1 = 0; % Default error structure: Normal distribution.
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
numstartpoints = 10; % Number of initial guesses for global optimization (Multistart).
B = 100;             % Number of bootstrap realizations for parameter uncertainty characterization.

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

flag1 = GLM;         % Selected growth model: Generalized Logistic Growth Model (GLM).
model_name1 = 'GLM'; % Name of the selected model.
fixI0 = 1;           % Boolean: Fix initial value to the first data point (true) or estimate it (false).

% <============================================================================>
% <=========== Parameters for Rolling Window Analysis =======================>
% <============================================================================>
% Settings for rolling window analysis:
% - windowsize1: Size of the moving window.
% - tstart1: Time point where rolling window analysis starts.
% - tend1: Time point where rolling window analysis ends.

windowsize1 = 20; % Size of the rolling window (e.g., 20 days).
tstart1 = 1;     % Start time point for rolling window analysis.
tend1 = 1;       % End time point for rolling window analysis.

end
