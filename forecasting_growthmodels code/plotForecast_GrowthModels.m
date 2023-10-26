function [param_doubling,seq_doublingtimes]=plotForecast_GrowthModels(tstart1_pass,tend1_pass,windowsize1_pass,forecastingperiod_pass)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% Fitting model to epidemic data with quantified uncertainty

close all

% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1 % Parameter estimation method


% <============================================================================>
% <=================== Load parameter values supplied by user =================>
% <============================================================================>

[cadfilename1_INP,caddisease_INP,datatype_INP, dist1_INP, numstartpoints_INP,M_INP,flag1_INP,model_name1_INP,fixI0_INP,getperformance_INP2,forecastingperiod_INP,windowsize1_INP,tstart1_INP,tend1_INP]=options_forecast;


% <============================================================================>
% <================================ Datasets properties ==============================>
% <============================================================================>

cadfilename1=cadfilename1_INP;

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

if strcmp('CUMULATIVE',upper(cadfilename1(1:10)))==1

    data(:,2)=[data(1,2);diff(data(:,2))]; % Incidence curve

end

% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP2; % flag or indicator variable (1/0) to calculate forecasting performance or not

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

% <==================================================================================>
% ============================ Rolling window analysis=====================================>
% <==================================================================================>

param_rs2=[];
param_as2=[];
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


factors=factor(length(tstart1:1:tend1));

if length(factors)==1
    rows=factors;
    cols=1;

elseif length(factors)==3
    rows=factors(1)*factors(2);
    cols=factors(3);
else
    rows=factors(1);
    cols=factors(2);
end


cc1=1;

for i=tstart1:1:tend1  %rolling window analysis

    load(strcat('./output/Forecast-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-',num2str(forecastingperiod),'.mat'))

    getperformance=getperformance_INP2; % flag or indicator variable (1/0) to calculate forecasting performance or not


    % <======================================================================================>
    % <======================= Plot parameter distributions and model fit and forecast ========================>
    % <======================================================================================>

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
    % <================================ Store parameter estimates =========================================>
    % <========================================================================================>

    param_rs2=[param_rs2; param_r];
    param_as2=[param_as2; param_a];
    param_ps2=[param_ps2; param_p];
    param_Ks2=[param_Ks2; param_K];
    param_I0s2=[param_I0s2; param_I0];
    param_alphas2=[param_alphas2; param_alpha];
    param_ds2=[param_ds2; param_d];

    % <========================================================================================>
    % <================================ Plot model fit and forecast ======================================>
    % <========================================================================================>

    subplot(2,4,[5 6 7 8])

    plot(timevect2,forecast_model12,'c')
    hold on

    % plot 95% PI

    LB1=quantile(forecast_model12',0.025)';
    LB1=(LB1>=0).*LB1;

    UB1=quantile(forecast_model12',0.975)';
    UB1=(UB1>=0).*UB1;

    median1=median(forecast_model12,2);

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

    % compute doubling times
    meandoublingtime=[];

    doublingtimess=zeros(30,M)+NaN;

    maxd=1;

    for j=1:M

        [tds,C0data,curve,doublingtimes]=getDoublingTimeCurve(fit_model1(:,j),DT,0);

        doublingtimess(1:length(doublingtimes),j)=doublingtimes;

        if maxd<length(doublingtimes)
            maxd=length(doublingtimes);
        end

        meandoublingtime=[meandoublingtime;mean(doublingtimes)];
    end

    doublingtimess=doublingtimess(1:maxd,1:M);

    seq_doublingtimes=[];

    for j=1:maxd

        index1=find(~isnan(doublingtimess(j,:)));

        seq_doublingtimes=[seq_doublingtimes;[j mean(doublingtimess(j,index1)) quantile(doublingtimess(j,index1),0.025) quantile(doublingtimess(j,index1),0.975) length(index1)./M]];

    end

    seq_doublingtimes % [ith doubling, mean, 95%CI LB, 95%CI UB, prob. i_th doubling]

    % Mean doubling times
    dmean=mean(meandoublingtime);
    dLB=quantile(meandoublingtime,0.025);
    dUB=quantile(meandoublingtime,0.975);

    param_doubling=[dmean dLB dUB]

    % <=============================================================================================>
    % <============================== Save file with doubling time estimates =======================>
    % <=============================================================================================>

    T = array2table(seq_doublingtimes);
    T.Properties.VariableNames(1:5) = {'i_th doubling','db mean','db 95%CI LB','db 95% CI UB','prob. i_th doubling'};
    writetable(T,strcat('./output/doublingtimes-flag1-',num2str(flag1),'-tstart-',num2str(i),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))



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


    %

    if length(tstart1:1:tend1)>1

        figure(400)

        subplot(rows,cols,cc1)

        plot(timevect2,forecast_model12,'c')
        hold on

        % plot 95% PI

        LB1=quantile(forecast_model12',0.025)';
        LB1=(LB1>=0).*LB1;

        UB1=quantile(forecast_model12',0.975)';
        UB1=(UB1>=0).*UB1;

        median1=median(forecast_model12,2);

        line1=plot(timevect2,median1,'r-')
        set(line1,'LineWidth',2)

        hold on
        line1=plot(timevect2,LB1,'r--')
        set(line1,'LineWidth',2)

        line1=plot(timevect2,UB1,'r--')
        set(line1,'LineWidth',2)

        % plot median model fit

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

        %xlabel('Time')
        %ylabel(strcat(caddisease,{' '},datatype))

        set(gca,'FontSize',16)
        set(gcf,'color','white')

        cc1=cc1+1;
    end

    %


    if getperformance & forecastingperiod>0 & length(data_all)<windowsize1+forecastingperiod

        [length(data_all) windowsize1+forecastingperiod]
        'entro'
        warning('Length of time series data is too short to evaluate the forecasting period indicated in <forecastingperiod>. Consider setting <getperformance> to 0 in options_forecast.m or extending the length of the time series.')

        forecastdata=[[timevect_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] [data_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] median1 LB1 UB1];

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
        writetable(T,strcat('./output/Forecast-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

    else

        if length(data_all)>=windowsize1+forecastingperiod
            forecastdata=[timevect_all(1:length(timevect1)+forecastingperiod) data_all(1:length(timevect1)+forecastingperiod) median1 LB1 UB1];
        else
            forecastdata=[[timevect_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] [data_all(1:end);zeros(windowsize1+forecastingperiod-length(data_all),1)+NaN] median1 LB1 UB1];
        end

        T = array2table(forecastdata);
        T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
        writetable(T,strcat('./output/Forecast-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

    end



    % <========================================================================================>
    % <========================================================================================>
    %                                  Plots forecasting performance metrics over predicted horizon
    % <========================================================================================>
    % <========================================================================================>

    if 1

        if forecastingperiod>0 & getperformance

            figure(200+i)

            subplot(2,2,1)

            line1=plot(MAEFS_model1(:,1),MAEFS_model1(:,2),'k-o')
            set(line1,'LineWidth',4)
            hold on

            xlabel('Forecasting horizon')
            ylabel('MAE')
            hold on
            set(gca,'FontSize', 24);
            set(gcf,'color','white')

            subplot(2,2,2)

            line1=plot(MSEFS_model1(:,1),MSEFS_model1(:,2),'k-o')
            set(line1,'LineWidth',4)
            hold on

            xlabel('Forecasting horizon')
            ylabel('MSE')
            hold on
            set(gca,'FontSize', 24);
            set(gcf,'color','white')

            subplot(2,2,3)

            line1=plot(PIFS_model1(:,1),PIFS_model1(:,2),'k-o')
            set(line1,'LineWidth',4)
            hold on

            xlabel('Forecasting horizon')
            ylabel('Coverage rate of the 95% PI')
            hold on
            set(gca,'FontSize', 24);
            set(gcf,'color','white')


            subplot(2,2,4)

            line1=plot(WISFS_model1(:,1),WISFS_model1(:,2),'k-o')
            set(line1,'LineWidth',4)
            hold on

            xlabel('Forecasting horizon')
            ylabel('WIS')
            hold on
            set(gca,'FontSize', 24);
            set(gcf,'color','white')
        end

        % <==================================================================================================>
        % <================================ Store performance metrics ==============================================>
        % <==================================================================================================>

        % store metrics for calibration
        RMSECSS2=[RMSECSS2;[RMSECS_model1(end,end)]];
        MSECSS2=[MSECSS2;[MSECS_model1(end,end)]];
        MAECSS2=[MAECSS2;[MAECS_model1(end,end)]];
        PICSS2=[PICSS2;[PICS_model1(end,end)]];
        MISCSS2=[MISCSS2;[MISCS_model1(end,end)]];

        WISCSS2=[WISCSS2;[WISC_model1(end,end)]];

        % store metrics for short-term forecasts
        if forecastingperiod>0 & getperformance

            RMSEFSS2=[RMSEFSS2;[RMSEFS_model1(end,end)]];
            MSEFSS2=[MSEFSS2;[MSEFS_model1(end,end)]];
            MAEFSS2=[MAEFSS2;[MAEFS_model1(end,end)]];
            PIFSS2=[PIFSS2;[PIFS_model1(end,end)]];
            MISFSS2=[MISFSS2;[MISFS_model1(end,end)]];

            WISFSS2=[WISFSS2;[WISFS_model1(end,end)]];

        end

        % <==================================================================================================>
        % <====================== Store quantiles of the calibration and forecasting periods and store ============================>
        % <==================================================================================================>

        quantilescs2=[quantilescs2;quantilesc];

        quantilesfs2=[quantilesfs2;quantilesf];


    end

end  % Rolling window analysis

if length(tstart1:1:tend1)>1

    figure(400)
    for c=1:cols
        subplot(rows,cols,(rows-1)*cols+c)
        xlabel('Time')
    end

    subplot(rows,cols,1)
    ylabel(strcat(caddisease,{' '},datatype))
    
end


% <==================================================================================================>
% <====================== Plot temporal variation of parameters from rolling window analysis ============================>
% <==================================================================================================>

if tend1>tstart1
    figure

    subplot(2,4,[1 2 3 4])

    plot(data(:,1),data(:,2),'ro-')
    xlabel('Time')
    ylabel('Cases')
    set(gca,'FontSize',24)
    set(gcf,'color','white')


    subplot(2,4,5)

    plot(tstart1:1:tend1,param_rs2(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_rs2(:,2),'b--')
    plot(tstart1:1:tend1,param_rs2(:,3),'b--')

    line1=plot(tstart1:1:tend1,smooth(param_rs2(:,1),5),'k--')
    set(line1,'LineWidth',3)


    ylabel('r')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

    subplot(2,4,6)
    plot(tstart1:1:tend1,param_as2(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_as2(:,2),'b--')
    plot(tstart1:1:tend1,param_as(:,3),'b--')

    line1=plot(tstart1:1:tend1,smooth(param_as2(:,1),5),'k--')
    set(line1,'LineWidth',3)

    ylabel('a')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

    subplot(2,4,7)
    plot(tstart1:1:tend1,param_ps2(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_ps2(:,2),'b--')
    plot(tstart1:1:tend1,param_ps2(:,3),'b--')

    line1=plot(tstart1:1:tend1,smooth(param_ps2(:,1),5),'k--')
    set(line1,'LineWidth',3)


    ylabel('p')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')

    subplot(2,4,8)
    plot(tstart1:1:tend1,param_Ks2(:,1),'ro-')
    hold on
    plot(tstart1:1:tend1,param_Ks2(:,2),'b--')
    plot(tstart1:1:tend1,param_Ks2(:,3),'b--')

    line1=plot(tstart1:1:tend1,smooth(param_Ks2(:,1),5),'k--')
    set(line1,'LineWidth',3)

    ylabel('K')
    set(gca,'FontSize',24)
    set(gcf,'color','white')
    xlabel('Time')
end

% <=====================================================================================================>
% <============================== Save file with AIC metrics ===========================================>
% <=====================================================================================================>

%[i AICc part1 part2 numparams]];

T = array2table(AICcs);
T.Properties.VariableNames(1:5) = {'time','AICc','AICc part1','AICc part2','numparams'};
writetable(T,strcat('./output/AICcs-rollingwindow-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))


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

performanceC=[(tstart1:1:tend1)' zeros(length(MAECSS2(:,1)),1)+windowsize1 MAECSS2(:,1)  MSECSS2(:,1) PICSS2(:,1) WISCSS2(:,1)];

T = array2table(performanceC);
T.Properties.VariableNames(1:6) = {'time','calibration_period','MAE','MSE','Coverage 95%PI','WIS'};
writetable(T,strcat('./output/performance-calibration-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

% <========================================================================================>
% <========================================================================================>
% <========================== Save csv file with forecasting performance metrics ============================>
% <========================================================================================>
% <========================================================================================>

if getperformance && forecastingperiod>0

    performanceF=[(tstart1:1:tend1)' zeros(length(MAEFSS2(:,1)),1)+forecastingperiod MAEFSS2(:,1)  MSEFSS2(:,1) PIFSS2(:,1) WISFSS2(:,1)];

    T = array2table(performanceF);
    T.Properties.VariableNames(1:6) = {'time','Horizon','MAE','MSE','Coverage 95%PI','WIS'};
    writetable(T,strcat('./output/performance-forecasting-flag1-',num2str(flag1),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(tstart1),'-tend-',num2str(tend1),'-calibrationperiod-',num2str(windowsize1),'-horizon-',num2str(forecastingperiod),'-',caddisease,'-',datatype,'.csv'))

end
