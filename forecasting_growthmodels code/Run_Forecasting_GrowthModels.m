function   AICcs=Run_Fit_GrowthModels(tstart1_pass,tend1_pass,windowsize1_pass,forecastingperiod_pass)

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

[cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP,flag1_INP,model_name1_INP,fixI0_INP,getperformance_INP,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP]=options_forecast;
   
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

flag1=flag1_INP; % Select the growth model

model_name1=model_name1_INP;

% <==================================================================================>
% <=============================== other parameters=======================================>
% <==================================================================================>

fixI0=fixI0_INP; % 0=Estimate the initial number of cases; 1 = Fix the initial number of cases according to the first data point in the time series


% <==============================================================================>
% <======================== Load epiemic data ========================================>
% <==============================================================================>

data=load(strcat('./input/',cadfilename1,'.txt'));

if strcmp('CUMULATIVE',upper(cadfilename1(1:10)))==1

    data(:,2)=[data(1,2);diff(data(:,2))]; % Incidence curve

end

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

if exist('forecastingperiod_pass','var')==1 & isempty(forecastingperiod_pass)==0

  forecastingperiod=forecastingperiod_pass; %forecast horizon (number of data points ahead)

else
    forecastingperiod=forecastingperiod_INP;

end

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

if exist('tstart1_pass','var')==1 & isempty(tstart1_pass)==0

   tstart1=tstart1_pass;

else
    tstart1=tstart1_INP;

end

if exist('tend1_pass','var')==1 & isempty(tend1_pass)==0

    tend1=tend1_pass;
else
    tend1=tend1_INP;

end

if exist('windowsize1_pass','var')==1 & isempty(windowsize1_pass)==0

    windowsize1=windowsize1_pass;
else
    windowsize1=windowsize1_INP;
end


% <===========================================================================================================>
% <====== Check that the number of estimated parameters is smaller than the number of data points= ===========>
% <===========================================================================================================>

numparams=get_nparams(method1,dist1,flag1,fixI0);

numparams
windowsize1

if numparams>=windowsize1

    error("Number of estimated parameters should be smaller than the calibration period. Consider increasing the length of the calibration period.")

end


% <==================================================================================>
% ============================ Rolling window analysis=====================================>
% <==================================================================================>

param_rs=[];
param_as=[];
param_ps=[];
param_Ks=[];
param_I0s=[];
param_alphas=[];
param_ds=[];

MCSES=[];

RMSECSS=[];
MSECSS=[];
MAECSS=[];
PICSS=[];
MISCSS=[];
RMSEFSS=[];
MSEFSS=[];
MAEFSS=[];
PIFSS=[];
MISFSS=[];

WISCSS=[];
WISFSS=[];

quantilescs=[];
quantilesfs=[];

if (tend1+windowsize1-1) > length(data(:,1))

    tend1= length(data(:,1))-windowsize1+1;
    'adjusting tend1'

end

AICcs=[];

for i=tstart1:1:tend1  %rolling window analysis

    figure(100+i)

    t_window=i:1:i+windowsize1-1;

    %calibration data
    datac=data(t_window,2);
    timevect1=DT*data(t_window,1);

    % calibration and future data
    data_all=data(i:end,2);
    timevect_all=DT*data(i:end,1);

    % first data point cannot be zero
    if datac(1)==0
        warning('first data point in the time series is zero')
        %continue
    end

    if getperformance & length(data_all)<windowsize1+forecastingperiod

        [length(data_all) windowsize1+forecastingperiod]
        
        error('Length of time series data is too short to evaluate the forecasting period indicated in <forecastingperiod>.')
        
    end


    I0=datac(1);

    data1=[data(t_window,1) datac];

    params0=initialParams(datac,flag1);


    [P_model1d,residual_model1, fitcurve_model1d, forecastcurve_model1, timevect2, initialguess,fval]=fit_model(data1,params0,1,numstartpoints,DT,flag1,0);


    %plot(timevect1,data1(:,2),'ko')
    %hold on
    %plot(timevect1, fitcurve_model1d,'r--')

    %xlabel('Time')
    %ylabel('Cases')

    [AICc,part1,part2,numparams]=getAICc(method1,dist1,flag1,fixI0,fval,length(data1(:,1)))

    AICcs=[AICcs;[i AICc part1 part2 numparams]];

    if (method1==0 & dist1==1) 

        factor1=1;

    elseif (method1==0 & dist1==2)  % estimate <factor1> which is the variance to mean ratio to scale the NB error structure

        binsize1=7;

        [ratios,~]=getMeanVarianceRatio(data1,binsize1,2);

        index1=find(ratios(:,1)>0);

        factor1=mean(ratios(index1,1));

        if factor1<1
            dist1=1;
            factor1=1;
        end

    elseif method1==0 & dist1==0 % normal distribution of the error structure

        var1=sum((fitcurve_model1d-datac).^2)./(length(fitcurve_model1d)-numparams); % last revised: 01 June 2022
        factor1=sqrt(var1);

    elseif method1==1
        dist1=1;
        factor1=1;

    elseif method1==3

        dist1=3;
        alpha=P_model1d(6); % VAR=mean+alpha*mean;
        factor1=alpha;

    elseif method1==4

        dist1=4;
        alpha=P_model1d(6); % VAR=mean+alpha*mean^2;
        factor1=alpha;

    elseif method1==5

        dist1=5;

        alpha=P_model1d(6); % VAR=mean+alpha*mean^d;
        d=P_model1d(7);

        factor1=alpha;

    end


    % <===========================================================================================>
    % <======= Derive parameter uncertainty of the best fitting models and save results ================================>
    % <===========================================================================================>


    fit_model1=[];
    forecast_model1=[];
    forecast_model12=[];

    Phatss_model1=zeros(M,7);

    f_model1_sims=[];

    for j=1:M


        f_model1_sim=AddErrorStructure(cumsum(fitcurve_model1d),1,dist1,factor1,d);

        f_model1_sims=[f_model1_sims f_model1_sim];


        % Fit model1 to bootstrap data
        data1=[timevect1 f_model1_sim];

        %params0=initialParams(data1(:,2),flag1);
        params0=P_model1d;

        [P_model1,residual_model1 fitcurve_model1 forecastcurve_model1 timevect2]=fit_model(data1,params0,fixI0,2,DT,flag1,forecastingperiod);

        fit_model1=[fit_model1 fitcurve_model1];

        forecast_model1=[forecast_model1 forecastcurve_model1];

        if method1==0 & dist1==0

            forecast_model12=[forecast_model12 AddErrorStructure(cumsum(forecastcurve_model1),20,dist1,factor1,0)];

        else

            forecast_model12=[forecast_model12 AddErrorStructure(cumsum(forecastcurve_model1),20,dist1,P_model1(end-1),P_model1(end))];

        end

        Phatss_model1(j,:)=P_model1;

    end %end bootstrapping loop

    data1=[timevect1 data1(:,2)];

    datalatest=data(i:end,1:2);

    % <==================================================================================================>
    % <========== Get forecast performance metrics for the model (if getperformance=1) =====================================>
    % <==================================================================================================>

    if 1

        [RMSECS_model1 MSECS_model1 MAECS_model1  PICS_model1 MISCS_model1 RMSEFS_model1 MSEFS_model1 MAEFS_model1 PIFS_model1 MISFS_model1]=computeforecastperformance(data1,datalatest,forecast_model1,forecast_model12,forecastingperiod);

        [WISC_model1,WISFS_model1]=computeWIS(data1,datalatest,forecast_model12,forecastingperiod)

        % store metrics for calibration
        RMSECSS=[RMSECSS;[RMSECS_model1(end,end)]];
        MSECSS=[MSECSS;[MSECS_model1(end,end)]];
        MAECSS=[MAECSS;[MAECS_model1(end,end)]];
        PICSS=[PICSS;[PICS_model1(end,end)]];
        MISCSS=[MISCSS;[MISCS_model1(end,end)]];

        WISCSS=[WISCSS;[WISC_model1(end,end)]];

        % store metrics for short-term forecasts
        if forecastingperiod>0 && isempty(RMSEFS_model1)==0

            RMSEFSS=[RMSEFSS;[RMSEFS_model1(end,end)]];
            MSEFSS=[MSEFSS;[MSEFS_model1(end,end)]];
            MAEFSS=[MAEFSS;[MAEFS_model1(end,end)]];
            PIFSS=[PIFSS;[PIFS_model1(end,end)]];
            MISFSS=[MISFSS;[MISFS_model1(end,end)]];

            WISFSS=[WISFSS;[WISFS_model1(end,end)]];

        end

    end

    % <==================================================================================================>
    % <========== Compute quantiles of the calibration and forecasting periods and store ================>
    % <==================================================================================================>

    [quantilesc,quantilesf]=computeQuantiles(data1,forecast_model12,forecastingperiod);

    quantilescs=[quantilescs;quantilesc];

    quantilesfs=[quantilesfs;quantilesf];


    % <========================================================================================>
    % <================================ Parameter estimates =========================================>
    % <========================================================================================>


    % estimate mean and 95% CI from distribution of parameter estimates
    param_r=[median(Phatss_model1(:,1)) quantile(Phatss_model1(:,1),0.025) quantile(Phatss_model1(:,1),0.975)];
    param_p=[median(Phatss_model1(:,2)) quantile(Phatss_model1(:,2),0.025) quantile(Phatss_model1(:,2),0.975)];

    param_a=[median(Phatss_model1(:,3)) quantile(Phatss_model1(:,3),0.025) quantile(Phatss_model1(:,3),0.975)];
    param_K=[median(Phatss_model1(:,4)) quantile(Phatss_model1(:,4),0.025) quantile(Phatss_model1(:,4),0.975)];

    param_I0=[median(Phatss_model1(:,5)) quantile(Phatss_model1(:,5),0.025) quantile(Phatss_model1(:,5),0.975)];

    param_alpha=[median(Phatss_model1(:,6)) quantile(Phatss_model1(:,6),0.025) quantile(Phatss_model1(:,6),0.975)];

    param_d=[median(Phatss_model1(:,7)) quantile(Phatss_model1(:,7),0.025) quantile(Phatss_model1(:,7),0.975)];

    MCSE=[std(Phatss_model1(:,1))/sqrt(M) std(Phatss_model1(:,2))/sqrt(M) std(Phatss_model1(:,3))/sqrt(M) std(Phatss_model1(:,4))/sqrt(M) std(Phatss_model1(:,5))/sqrt(M) std(Phatss_model1(:,6))/sqrt(M) std(Phatss_model1(:,7))/sqrt(M)];

    param_rs=[param_rs; param_r];
    param_as=[param_as; param_a];
    param_ps=[param_ps; param_p];
    param_Ks=[param_Ks; param_K];
    param_I0s=[param_I0s; param_I0];
    param_alphas=[param_alphas; param_alpha];
    param_ds=[param_ds; param_d];


    MCSES=[MCSES;MCSE];

    cad1=strcat('r=',num2str(param_r(end,1),2),' (95% CI:',num2str(param_r(end,2),2),',',num2str(param_r(end,3),2),')')
    cad2=strcat('p=',num2str(param_p(end,1),2),' (95% CI:',num2str(param_p(end,2),2),',',num2str(param_p(end,3),2),')')

    cad3=strcat('a=',num2str(param_a(end,1),2),' (95% CI:',num2str(param_a(end,2),2),',',num2str(param_a(end,3),2),')')
    cad4=strcat('K=',num2str(param_K(end,1),2),' (95% CI:',num2str(param_K(end,2),2),',',num2str(param_K(end,3),2),')')

    cad5=strcat('alpha=',num2str(param_alpha(end,1),2),' (95% CI:',num2str(param_alpha(end,2),2),',',num2str(param_alpha(end,3),2),')')

    cad6=strcat('d=',num2str(param_d(end,1),2),' (95% CI:',num2str(param_d(end,2),2),',',num2str(param_d(end,3),2),')')


    % Plot results

    LB1=quantile(forecast_model12',0.025)';
    LB1=(LB1>=0).*LB1;

    UB1=quantile(forecast_model12',0.975)';
    UB1=(UB1>=0).*UB1;

    median1=median(forecast_model12,2);

    if 1


        % <========================================================================================>
        % <======================= Plot empirical distributions of the parameters ========================>
        % <========================================================================================>

        figure(100+i)
        subplot(2,4,1)
        hist(Phatss_model1(:,1))
        hold on

        line2=[param_r(1,2) 10;param_r(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('r')
        ylabel('Frequency')

        title(cad1)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,4,2)
        hist(Phatss_model1(:,2))
        hold on

        line2=[param_p(1,2) 10;param_p(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('p')
        ylabel('Frequency')

        title(cad2)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,4,3)

        hist(Phatss_model1(:,3))
        hold on

        line2=[param_a(1,2) 10;param_a(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('a')
        ylabel('Frequency')
        title(cad3)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,4,4)
        hist(Phatss_model1(:,4))
        hold on

        line2=[param_K(1,2) 10;param_K(1,3) 10];
        line1=plot(line2(:,1),line2(:,2),'r--')
        set(line1,'LineWidth',2)

        xlabel('K')
        ylabel('Frequency')

        title(cad4)

        set(gca,'FontSize', 24);
        set(gcf,'color','white')


        % <========================================================================================>
        % <================================ Plot model fit and forecast ======================================>
        % <========================================================================================>

        subplot(2,4,[5 6 7 8])

        plot(timevect2,forecast_model12,'c')
        hold on

        % plot 95% PI


        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)

        hold on
        line1=plot(timevect2,LB1,'r--')
        set(line1,'LineWidth',2)

        line1=plot(timevect2,UB1,'r--')
        set(line1,'LineWidth',2)

        % plot mean model fit

        color1=gray(8);
        line1=plot(timevect1,fit_model1,'color',color1(6,:))
        set(line1,'LineWidth',1)

        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)

        % plot the data

        line1=plot(timevect_all,data_all,'bo')
        set(line1,'LineWidth',2)

        line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast_model12',0.975))*1.5];

        if forecastingperiod>0
            line1=plot(line2(:,1),line2(:,2),'k--')
            set(line1,'LineWidth',2)
        end

        axis([timevect1(1) timevect2(end) 0 max(quantile(forecast_model12',0.975))*1.5])

        xlabel('Time')
        ylabel(strcat(caddisease,{' '},datatype))

        set(gca,'FontSize',24)
        set(gcf,'color','white')

        title(model_name1)
    end


    % <=========================================================================================>
    % <================================ Save short-term forecast results ==================================>
    % <=========================================================================================>

    save(strcat('./output/Forecast-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'-mat')

end % rolling window analysis

%save model parameters from tstart1 to tend1
save(strcat('./output/parameters-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'), 'param_rs', 'param_as', 'param_ps', 'param_Ks', 'param_I0s', 'param_alphas', 'param_ds','-mat')

%save calibration performance metrics from tstart1 to tend1 (AICc, MSE,
%MAE, coverage 95% PI, MIS, WIS)
save(strcat('./output/performanceCalibration-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'AICcs','RMSECSS', 'MSECSS', 'MAECSS', 'PICSS','MISCSS','WISCSS','-mat')

%save quantiles of the model fits from tstart1 to tend1
save(strcat('./output/QuantilesCalibration-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'quantilescs','-mat')

% <=====================================================================================================>
% <============================== Save file with AIC metrics ===========================================>
% <=====================================================================================================>

%[i AICc part1 part2 numparams]];

T = array2table(AICcs);
T.Properties.VariableNames(1:5) = {'time','AICc','AICc part1','AICc part2','numparams'};
writetable(T,strcat('./output/AICcs-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


if forecastingperiod>0

    %save forecasting performance metrics from tstart1 to tend1 (AICc, MSE,
    %MAE, coverage 95% PI, MIS, WIS)
    save(strcat('./output/performanceForecasting-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'RMSEFSS', 'MSEFSS', 'MAEFSS', 'PIFSS','MISFSS','WISFSS','-mat')

    %save quantiles of the model forecasts from tstart1 to tend1
    save(strcat('./output/QuantilesForecastingPerformance-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'),'quantilesfs','-mat')

end


% <=============================================================================================>
% <================= Save csv file with parameters from rolling window analysis ====================================>
% <=============================================================================================>

if method1==3 | method1==4  %save parameter alpha. VAR=mean+alpha*mean; VAR=mean+alpha*mean^2;

    rollparams=[(tstart1:1:tend1)' param_rs(:,1:end) param_ps(:,1:end) param_as(:,1:end) param_Ks(:,1:end) param_I0s(:,1:end) param_alphas(:,1:end)];
    T = array2table(rollparams);
    T.Properties.VariableNames(1:19) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','a mean','a LB','a UB','K0 mean','K0 LB','K0 UB','I0 mean','I0 LB','I0 UB','alpha mean','alpha LB','alpha UB'};

    rollparams=[(tstart1:1:tend1)' MCSES(:,1:6)];
    T2 = array2table(rollparams);
    T2.Properties.VariableNames(1:7) = {'time','r MCSE','p MCSE','a MCSE','K0 MCSE','I0 MCSE','alpha MCSE'};

elseif method1==5

    rollparams=[(tstart1:1:tend1)' param_rs(:,1:end) param_ps(:,1:end) param_as(:,1:end) param_Ks(:,1:end) param_I0s(:,1:end) param_alphas(:,1:end) param_ds(:,1:end)];
    T = array2table(rollparams);
    T.Properties.VariableNames(1:22) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','a mean','a LB','a UB','K0 mean','K0 LB','K0 UB','I0 mean','I0 LB','I0 UB','alpha mean','alpha LB','alpha UB','d mean','d LB','d UB'};

    rollparams=[(tstart1:1:tend1)' MCSES(:,1:7)];
    T2 = array2table(rollparams);
    T2.Properties.VariableNames(1:8) = {'time','r MCSE','p MCSE','a MCSE','K0 MCSE','I0 MCSE','alpha MCSE','d MCSE'};

else

    rollparams=[(tstart1:1:tend1)' param_rs(:,1:end) param_ps(:,1:end) param_as(:,1:end) param_Ks(:,1:end) param_I0s(:,1:end)];
    T = array2table(rollparams);
    T.Properties.VariableNames(1:16) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','a mean','a LB','a UB','K0 mean','K0 LB','K0 UB','I0 mean','I0 LB','I0 UB'};

    rollparams=[(tstart1:1:tend1)' MCSES(:,1:5)];
    T2 = array2table(rollparams);
    T2.Properties.VariableNames(1:6) = {'time','r MCSE','p MCSE','a MCSE','K0 MCSE','I0 MCSE'};

end

writetable(T,strcat('./output/parameters-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

writetable(T2,strcat('./output/MCSES-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

% <========================================================================================>
% <========================================================================================>
% <========================== Save csv file with calibration performance metrics ============================>
% <========================================================================================>
% <========================================================================================>

performanceC=[(tstart1:1:tend1)' zeros(length(MAECSS(:,1)),1)+windowsize1 MAECSS(:,1)  MSECSS(:,1) PICSS(:,1) WISCSS(:,1)];

T = array2table(performanceC);
T.Properties.VariableNames(1:6) = {'time','calibration_period','MAE','MSE','Coverage 95%PI','WIS'};
writetable(T,strcat('./output/performance-calibration-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

% <========================================================================================>
% <========================================================================================>
% <========================== Save csv file with forecasting performance metrics ============================>
% <========================================================================================>
% <========================================================================================>

if getperformance && forecastingperiod>0 && isempty(length(MAEFSS))==0

    performanceF=[(tstart1:1:tend1)' zeros(length(MAEFSS(:,1)),1)+forecastingperiod MAEFSS(:,1)  MSEFSS(:,1) PIFSS(:,1) WISFSS(:,1)];

    T = array2table(performanceF);
    T.Properties.VariableNames(1:6) = {'time','Horizon','MAE','MSE','Coverage 95%PI','WIS'};
    writetable(T,strcat('./output/performance-forecasting-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

end


