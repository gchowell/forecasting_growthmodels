function [RMSECS MSECS MAECS  PICS MISCS RMSEFS MSEFS MAEFS PIFS MISFS]=computeforecastperformance(data1,datalatest,curvesforecasts1,curvesforecasts2,forecastingperiod)

RMSECS=[];
MSECS=[];
MAECS=[];
PICS=[];
MISCS=[];

RMSEFS=[];
MSEFS=[];
MAEFS=[];
PIFS=[];
MISFS=[];

% **** compute performance metrics ******

% calibration period

tf2=length(data1(:,1));

dataCperiod=1:tf2;

%ycmean=plims(curvesforecasts1(1:tf2,:)',0.5)';
ycmean=median(curvesforecasts1(1:tf2,:),2);

datac=datalatest(dataCperiod,2);

RMSEc=sqrt(mean((datac-ycmean).^2));

MSEc=mean((datac-ycmean).^2);
MAEc=mean(abs(datac-ycmean));

% 95% prediction coverage
Lt=quantile(curvesforecasts2(1:tf2,:)',0.025)';
Lt=(Lt>=0).*Lt;

Ut=quantile(curvesforecasts2(1:tf2,:)',0.975)';
Ut=(Ut>=0).*Ut;

coverage1=find(datac>=Lt & datac<=Ut);

'prediction coverage (%)='
PIc=100*(length(coverage1)./length(datac));

PICS=[PICS;[  PIc]];

MISc=mean([(Ut-Lt)+(2/0.05).*(Lt-datac).*(datac<Lt)+(2/0.05)*(datac-Ut).*(datac>Ut)]);


%figure
%plot(WIS)
%pause

%

MISCS=[MISCS;[  MISc]];

RMSECS=[RMSECS;[  RMSEc]];

MSECS=[MSECS;[  MSEc]];

MAECS=[MAECS;[  MAEc]];

% forecasting period %%%%%
RMSEFS=[]; MSEFS=[]; MAEFS=[]; PIFS=[]; MISFS=[];

if forecastingperiod>0
    
    for forecastingperiod1=1:1:forecastingperiod
        
        dataFperiod=tf2+1:tf2+forecastingperiod1;
        
        
        %yfmean=plims(curvesforecasts1(tf2+1:end,:)',0.5)';
        %yfmean=mean(curvesforecasts1(tf2+1:end,:),2);
        yfmean=median(curvesforecasts1(dataFperiod,:),2);
        
        %yf=curvesforecasts2(tf2+1:end,:);
        
        datagroundtruth=datalatest;
        
        dataf=datagroundtruth(dataFperiod,2);
        
        
        %plot ground truth data
        %line1=plot(datagroundtruth(1:end-1,1)*DT,datagroundtruth(1:end-1,2),'ko')
        %set(line1,'LineWidth',2)
        
        max1=max(datagroundtruth(:,2));
        
        %axis([timevect2(1) timevect2(end)+7*3 0 max1+max1*0.2])
        
        
        RMSEf=sqrt(mean((dataf-yfmean).^2));
        
        MSEf=mean((dataf-yfmean).^2);
        MAEf=mean(abs(dataf-yfmean));
        
        % 95% prediction coverage
        Lt=quantile(curvesforecasts2(dataFperiod,:)',0.025)';
        Lt=(Lt>=0).*Lt;

        Ut=quantile(curvesforecasts2(dataFperiod,:)',0.975)';
        Ut=(Ut>=0).*Ut;

        coverage1=find(dataf>=Lt & dataf<=Ut);
        
        PIf=100*(length(coverage1)./length(dataf));
        
        PIFS=[PIFS;[  forecastingperiod1 PIf]];
        
        %Mean Interval Score (MIS)
        
        MISf=mean([(Ut-Lt)+(2/0.05).*(Lt-dataf).*(dataf<Lt)+(2/0.05)*(dataf-Ut).*(dataf>Ut)]);
        
        MISFS=[MISFS;[  forecastingperiod1 MISf]];
        
        RMSEFS=[RMSEFS;[  forecastingperiod1 RMSEf]];
        
        MSEFS=[MSEFS;[  forecastingperiod1 MSEf]];
        
        MAEFS=[MAEFS;[  forecastingperiod1 MAEf]];
        
    end
    
end

