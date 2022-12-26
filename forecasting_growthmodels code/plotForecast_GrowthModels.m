% Fitting model to epidemic data with quantified uncertainty

clear
clear global
close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method


% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

[cadfilename1_INP,caddisease_INP,datatype_INP, DT_INP, dist1_INP, numstartpoints_INP,M_INP,flag1_INP,model_name1_INP,fixI0_INP,getperformance_INP,forecastingperiod_INP, printscreen1_INP,windowsize1_INP,tstart1_INP,tend1_INP]=options_forecast;


% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

DT=DT_INP;

caddisease=caddisease_INP;
datatype=datatype_INP;

% <=============================================================================>
% <=========================== Parameter estimation ============================>
% <=============================================================================>

method1=0; % Type of estimation method: 0 = LSQ

d=1;

dist1=dist1_INP; %Define dist1 which is the type of error structure:

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squard=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

% % Define dist1 which is the type of error structure:
% switch method1
%
%     case 0
%
%         dist1=0; % Normnal distribution to model error structure
%
%         %dist1=2; % error structure type (Poisson=1; NB=2)
%
%         %factor1=1; % scaling factor for VAR=factor1*mean
%
%
%     case 3
%         dist1=3; % VAR=mean+alpha*mean;
%
%     case 4
%         dist1=4; % VAR=mean+alpha*mean^2;
%
%     case 5
%         dist1=5; % VAR=mean+alpha*mean^d;
%
% end

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

flag1=flag1_INP; % Sequence of subepidemic growth models considered in epidemic trajectory

model_name1=model_name1_INP;

% <==================================================================================>
% <=============================== other parameters=======================================>
% <==================================================================================>

fixI0=fixI0_INP; % 0=Estimate the initial number of cases; 1 = Fix the initial number of cases according to the first data point in the time series


% <==============================================================================>
% <======================== Load epiemic data ========================================>
% <==============================================================================>

data=load(strcat('./input/',cadfilename1,'.txt'));


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

forecastingperiod=forecastingperiod_INP; %forecast horizon (number of data points ahead)

printscreen1=printscreen1_INP;  % print plots with the results

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=windowsize1_INP;  %moving window size
tstart1=tstart1_INP; % time of start of rolling window analysis
tend1=tend1_INP;  %time end of the rolling window analysis
%tend1=length(data(:,1));

% <==================================================================================>
% ============================ Rolling window analysis=====================================>
% <==================================================================================>

param_rs2=[];
param_ps2=[];
param_Ks2=[];
param_I0s2=[];
param_alphas2=[];
param_ds2=[];

RMSECSS2=[];
MSECSS2=[];
MAECSS2=[];
PICSS2=[];
MISCSS2=[];
RMSEFSS2=[];
MSEFSS2=[];
MAEFSS2=[];
PIFSS2=[];
MISFSS2=[];

WISCSS2=[];
WISFSS2=[];

quantilescs2=[];
quantilesfs2=[];

% <=========================================================================================>
% <================================ Load short-term forecast results ==================================>
% <=========================================================================================>

for i=tstart1:1:tend1-windowsize1+1  %rolling window analysis

    load(strcat('./output/Forecast-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'))

    % <======================================================================================>
    % <======================= Plot parameter distributions and model fit and forecast ========================>
    % <======================================================================================>

    figure(100+i)
    subplot(2,3,1)
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

    subplot(2,3,2)
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

    subplot(2,3,3)
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
    % <================================ Store parameter estimates =========================================>
    % <========================================================================================>

    param_rs2=[param_rs2; param_r];
    param_ps2=[param_ps2; param_p];
    param_Ks2=[param_Ks2; param_K];
    param_I0s2=[param_I0s2; param_I0];
    param_alphas2=[param_alphas2; param_alpha];
    param_ds2=[param_ds2; param_d];

    % <========================================================================================>
    % <================================ Plot model fit and forecast ======================================>
    % <========================================================================================>

    subplot(2,3,[4 5 6])

    plot(timevect2,forecast_model12,'c')
    hold on

    % plot 95% PI

    LB1=quantile(forecast_model12',0.025)';
    LB1=(LB1>=0).*LB1;

    UB1=quantile(forecast_model12',0.975)';
    UB1=(UB1>=0).*UB1;

    mean1=mean(forecast_model12,2);

    line1=plot(timevect2,mean1,'r-')
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

    % plot the data

    line1=plot(timevect_all,data_all,'bo')
    set(line1,'LineWidth',2)

    line2=[timevect1(end) 0;timevect1(end) max(quantile(forecast_model12',0.975))*1.5];

    if forecastingperiod>0
        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)
    end

    axis([timevect1(1) timevect2(end) 0 max(quantile(forecast_model12',0.975))*1.5])

    xlabel('Time (days)')
    ylabel(strcat(caddisease,{' '},datatype))

    set(gca,'FontSize',24)
    set(gcf,'color','white')

    title(model_name1)

    if getperformance & length(data_all)<windowsize1+forecastingperiod

        [length(data_all) windowsize1+forecastingperiod]
        'entro'
        warning('Length of time series data is too short to evaluate the forecasting period indicated in <forecastingperiod>. Consider setting <getperformance> to 0 in options_forecast.m or extending the length of the time series.')

        forecastdata=[[timevect_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] [data_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] median1 LB1 UB1];

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
        writetable(T,strcat('./output/Forecast-tstart-',num2str(tstart1),'-',caddisease,'-',datatype,'.csv'))

    else

        if length(data_all)>=windowsize1+forecastingperiod
            forecastdata=[timevect_all(1:length(timevect1)+forecastingperiod) data_all(1:length(timevect1)+forecastingperiod) median1 LB1 UB1];
        else
            forecastdata=[[timevect_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] [data_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] median1 LB1 UB1];
        end

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
        writetable(T,strcat('./output/Forecast-tstart-',num2str(tstart1),'-',caddisease,'-',datatype,'.csv'))

    end



    % <========================================================================================>
    % <========================================================================================>
    %                                  Plots forecasting performance metrics over predicted horizon
    % <========================================================================================>
    % <========================================================================================>

    if getperformance

        figure(200+i)

        subplot(2,2,1)

        line1=plot(MAEFS_model1(:,1),MAEFS_model1(:,2),'k')
        set(line1,'LineWidth',4)
        hold on

        xlabel('Forecasting horizon (days)')
        ylabel('MAE')
        hold on
        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,2,2)

        line1=plot(MSEFS_model1(:,1),MSEFS_model1(:,2),'k')
        set(line1,'LineWidth',4)
        hold on

        xlabel('Forecasting horizon (days)')
        ylabel('MSE')
        hold on
        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        subplot(2,2,3)

        line1=plot(PIFS_model1(:,1),PIFS_model1(:,2),'k')
        set(line1,'LineWidth',4)
        hold on

        xlabel('Forecasting horizon (days)')
        ylabel('Coverage rate of the 95% PI')
        hold on
        set(gca,'FontSize', 24);
        set(gcf,'color','white')


        subplot(2,2,4)

        line1=plot(WISFS_model1(:,1),WISFS_model1(:,2),'k')
        set(line1,'LineWidth',4)
        hold on

        xlabel('Forecasting horizon (days)')
        ylabel('WIS')
        hold on
        set(gca,'FontSize', 24);
        set(gcf,'color','white')

        % <========================================================================================>
        % <========================================================================================>
        % <========================== Save file with forecasting performance metrics ============================>
        % <========================================================================================>
        % <========================================================================================>

        performance=[MAEFS_model1(:,1:2)  MSEFS_model1(:,2) PIFS_model1(:,2) WISFS_model1(:,2)];

        T = array2table(performance);
        T.Properties.VariableNames(1:5) = {'Horizon','MAE','MSE','Coverage 95%PI','WIS'};
        writetable(T,strcat('./output/performance-forecasting-horizon-',caddisease,'-',datatype,'.csv'))

        % <==================================================================================================>
        % <================================ Store performance metrics ==============================================>
        % <==================================================================================================>

        % store metrics for calibration
        RMSECSS2=[RMSECSS2;[i RMSECS_model1(end,end)]];
        MSECSS2=[MSECSS2;[i MSECS_model1(end,end)]];
        MAECSS2=[MAECSS2;[i MAECS_model1(end,end)]];
        PICSS2=[PICSS2;[i PICS_model1(end,end)]];
        MISCSS2=[MISCSS2;[i MISCS_model1(end,end)]];

        WISCSS2=[WISCSS2;[i WISC_model1(end,end)]];

        % store metrics for short-term forecasts
        if forecastingperiod>0

            RMSEFSS2=[RMSEFSS2;[i RMSEFS_model1(end,end)]];
            MSEFSS2=[MSEFSS2;[i MSEFS_model1(end,end)]];
            MAEFSS2=[MAEFSS2;[i MAEFS_model1(end,end)]];
            PIFSS2=[PIFSS2;[i PIFS_model1(end,end)]];
            MISFSS2=[MISFSS2;[i MISFS_model1(end,end)]];

            WISFSS2=[WISFSS2;[i WISFS_model1(end,end)]];

        end

        % <==================================================================================================>
        % <====================== Store quantiles of the calibration and forecasting periods and store ============================>
        % <==================================================================================================>

        quantilescs2=[quantilescs2;quantilesc];

        quantilesfs2=[quantilesfs2;quantilesf];


    end

end  % Rolling window analysis

% <==================================================================================================>
% <====================== Plot temporal variation of parameters from rolling window analysis ============================>
% <==================================================================================================>

figure

subplot(1,3,1)
plot(tstart1:1:tend1-windowsize1+1,param_rs(:,1),'ro-')
hold on
plot(tstart1:1:tend1-windowsize1+1,param_rs(:,2),'b--')
plot(tstart1:1:tend1-windowsize1+1,param_rs(:,3),'b--')

ylabel('r')
set(gca,'FontSize',24)
set(gcf,'color','white')
xlabel('Time')

subplot(1,3,2)
plot(tstart1:1:tend1-windowsize1+1,param_ps(:,1),'ro-')
hold on
plot(tstart1:1:tend1-windowsize1+1,param_ps(:,2),'b--')
plot(tstart1:1:tend1-windowsize1+1,param_ps(:,3),'b--')

ylabel('p')
set(gca,'FontSize',24)
set(gcf,'color','white')
xlabel('Time')

subplot(1,3,3)
plot(tstart1:1:tend1-windowsize1+1,param_Ks(:,1),'ro-')
hold on
plot(tstart1:1:tend1-windowsize1+1,param_Ks(:,2),'b--')
plot(tstart1:1:tend1-windowsize1+1,param_Ks(:,3),'b--')

ylabel('K')
set(gca,'FontSize',24)
set(gcf,'color','white')
xlabel('Time')


