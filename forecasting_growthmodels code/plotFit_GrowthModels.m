% plot model fit to epidemic data with quantified uncertainty

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

[cadfilename1_INP,caddisease_INP,datatype_INP,DT_INP, dist1_INP, numstartpoints_INP,M_INP,flag1_INP,model_name1_INP,fixI0_INP,printscreen1_INP,windowsize1_INP,tstart1_INP,tend1_INP]=options_fit;


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

printscreen1=printscreen1_INP;  % print plots with the results

% <==================================================================================>
% <========================== Parameters of the rolling window analysis =========================>
% <==================================================================================>

windowsize1=windowsize1_INP;  %moving window size
tstart1=tstart1_INP; % time of start of rolling window analysis
tend1=tend1_INP;  %time end of the rolling window analysis
%tend1=length(data(:,1));


% <=========================================================================================>
% <================================ Load short-term forecast results ==================================>
% <=========================================================================================>

for i=tstart1:1:tend1-windowsize1+1  %rolling window analysis

    load(strcat('./output/Forecast-growthModel-',cadfilename1,'-flag1-',num2str(flag1(1)),'-fixI0-',num2str(fixI0),'-method-',num2str(method1),'-dist-',num2str(dist1),'-tstart-',num2str(i),'-calibrationperiod-',num2str(windowsize1),'-forecastingperiod-0.mat'))

    % <======================================================================================>
    % <======================= Plot parameter distributions and model fit and forecast ========================>
    % <======================================================================================>

    figure(101+i)
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

    median1=median(forecast_model12,2);

    line1=plot(timevect2,median1,'r-')
    set(line1,'LineWidth',2)

    hold on
    line1=plot(timevect2,LB1,'r--')
    set(line1,'LineWidth',2)

    line1=plot(timevect2,UB1,'r--')
    set(line1,'LineWidth',2)

    % plot model fit

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


    fitdata=[timevect1 data_all(1:length(timevect1)) median1 LB1 UB1];

    T = array2table(fitdata);
    T.Properties.VariableNames(1:5) = {'time','data','median','LB','UB'};
    writetable(T,strcat('./output/Fit-i-',num2str(i),'-',caddisease,'-',datatype,'.csv'))

    
end

% <==================================================================================================>
% <====================== Plot temporal variation of parameters from rolling window analysis ============================>
% <==================================================================================================>

figure

subplot(2,3,[1 2 3])

plot(data(:,1),data(:,2),'ro-')
xlabel('Time')
ylabel('Cases')
set(gca,'FontSize',24)
set(gcf,'color','white')


subplot(2,3,4)

plot(tstart1:1:tend1-windowsize1+1,param_rs(:,1),'ro-')
hold on
plot(tstart1:1:tend1-windowsize1+1,param_rs(:,2),'b--')
plot(tstart1:1:tend1-windowsize1+1,param_rs(:,3),'b--')

line1=plot(tstart1:1:tend1-windowsize1+1,smooth(param_rs(:,1),5),'k--')
set(line1,'LineWidth',3)


ylabel('r')
set(gca,'FontSize',24)
set(gcf,'color','white')
xlabel('Time')

subplot(2,3,5)
plot(tstart1:1:tend1-windowsize1+1,param_ps(:,1),'ro-')
hold on
plot(tstart1:1:tend1-windowsize1+1,param_ps(:,2),'b--')
plot(tstart1:1:tend1-windowsize1+1,param_ps(:,3),'b--')

line1=plot(tstart1:1:tend1-windowsize1+1,smooth(param_ps(:,1),5),'k--')
set(line1,'LineWidth',3)


ylabel('p')
set(gca,'FontSize',24)
set(gcf,'color','white')
xlabel('Time')

subplot(2,3,6)
plot(tstart1:1:tend1-windowsize1+1,param_Ks(:,1),'ro-')
hold on
plot(tstart1:1:tend1-windowsize1+1,param_Ks(:,2),'b--')
plot(tstart1:1:tend1-windowsize1+1,param_Ks(:,3),'b--')

line1=plot(tstart1:1:tend1-windowsize1+1,smooth(param_Ks(:,1),5),'k--')
set(line1,'LineWidth',3)


ylabel('K')
set(gca,'FontSize',24)
set(gcf,'color','white')
xlabel('Time')

% <=============================================================================================>
% <================= Save file with parameters from rolling window analysis ====================================>
% <=============================================================================================>

performance=[(tstart1:1:tend1-windowsize1+1)' param_rs(:,1:end) param_ps(:,1:end) param_Ks(:,1:end)];

T = array2table(performance);
T.Properties.VariableNames(1:10) = {'time','r mean','r LB','r UB','p mean','p LB','p UB','K0 mean','K0 LB','K0 UB'};
writetable(T,strcat('./output/parameters-rollingwindow-',caddisease,'-',datatype,'.csv'))

