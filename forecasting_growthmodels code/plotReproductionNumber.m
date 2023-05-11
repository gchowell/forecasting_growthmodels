% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>


% Plot model fits and effective reproduction number during the calibration
% and forecasting periods derived from the top-ranked models

clear
close all

% <============================================================================>
% <=================== Declare global variables =======================================>
% <============================================================================>

global invasions
global timeinvasions
global Cinvasions
global npatches_fixed
global onset_fixed

global method1 dist1 factor1

global smoothfactor1

global calibrationperiod1

% <============================================================================>
% <================== Load the parameter values ===============================>
% <============================================================================>

% options.m
[outbreakx_INP, caddate1_INP, cadregion_INP, caddisease_INP, datatype_INP, DT_INP, datevecfirst1_INP, datevecend1_INP, numstartpoints_INP, topmodelsx_INP, M_INP, flag1_INP]=options

% options_forecast.m
[getperformance_INP, deletetempfiles_INP, forecastingperiod_INP, printscreen1_INP, weight_type1_INP]=options_forecast

%options_Rt.m
[type_GId1_INP,mean_GI1_INP,var_GI1_INP]=options_Rt

% <============================================================================>
% <================================ Dataset ===================================>
% <============================================================================>

outbreakx=outbreakx_INP;

caddate1=caddate1_INP;

cadregion=cadregion_INP;

caddisease=caddisease_INP;

datatype=datatype_INP;

datevecfirst1=datevecfirst1_INP;

datevecend1=datevecend1_INP;


DT=DT_INP; % temporal resolution in days (1=daily data, 7=weekly data, 365=yearly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
elseif DT==365
    cadtemporal='yearly';
end


cadfilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1);


% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

%smoothfactor1=7; % <smoothfactor1>-day rolling average smoothing of the case series

%calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)

% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

%method1=0;  % Type of estimation method: 0 = LSQ

%dist1=0; % Normnal distribution to model error structure

M=M_INP; %number of bootstrap realizations to generate uncertainty estimates

% <==============================================================================>
% <========================= Growth model ==========================================>
% <==============================================================================>

npatchess2=npatches_fixed;  % maximum number of subepidemics considered in epidemic trajectory fit

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flagss2=flag1_INP; % Sequence of subepidemic growth models considered in epidemic trajectory


% <===============================================================================================>
% <============= Number of best fitting models used to generate ensemble model ===================>
% <===============================================================================================>

topmodels1=1:topmodelsx_INP;

if npatches_fixed==1
    topmodels1=1;
end

factors=factor(length(topmodels1));
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


% <==============================================================================>
% <========================== Forecasting parameters ===================================>
% <==============================================================================>

getperformance=getperformance_INP; % flag or indicator variable (1/0) to calculate forecasting performance or not

deletetempfiles=deletetempfiles_INP; %flag or indicator variable (1/0) to delete Forecast..mat files after use

forecastingperiod=forecastingperiod_INP; %forecast horizon (number of data points ahead)

printscreen1=printscreen1_INP;  % print plots with the results

% <==============================================================================>
% <====================== weighting scheme for ensemble model ============================>
% <==============================================================================>

weight_type1=weight_type1_INP; % -1= equally weighted from the top models, 0=based on AICc, 1= based on relative likelihood (Akaike weights), 2=based on WISC during calibration, 3=based on WISF during forecasting performance at previous time period (week)

WISC_hash=zeros(length(topmodels1),1); % vector that saves the WISC based on calibration to be used with weight_type1=2

WISF_hash=zeros(length(topmodels1),200); % vector that saves the WISF based on prior forecasting performance to be used with weight_type1=2


% <=======================================================================================>
% <========================== Reproduction number number parameters =======================>
% <========================================================================================>

type_GId1=type_GId1_INP; % type of GI distribution 1=Gamma, 2=Exponential, 3=Delta

mean_GI1=mean_GI1_INP;  % mean of the generation interval distribution

var_GI1=var_GI1_INP; % variance of the generation interval distribution




RMSES=[];
MAES=[];
MSES=[];
PIS=[];
MISS=[];

cc1=1;

AICc_bests=[];


param_Rt=[];

for run_id=-1

    cc1=1;

    for rankx=topmodels1

        % <========================================================================================>
        % <================================ Load model results ====================================>
        % <========================================================================================>

        load (strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-0-smoothing-',num2str(smoothfactor1),'-',cadfilename2,'-flag1-',num2str(flagss2(1)),'-flag1-',num2str(flagss2(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rankx),'.mat'))

        d_hat=1;

        cc3=1;


        % <========================================================================================>
        % <================================ Compute short-term forecast ===========================>
        % <========================================================================================>


        timevect=(data1(:,1));


        timevect2=(0:t_window(end)-1+forecastingperiod);

        % vector to store forecast mean curves
        curvesforecasts1=[];

        % vector to store forecast prediction curves
        curvesforecasts2=[];

        color1=['r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';'r-';'b-';'g-';'m-';'c-';'k-';'y-';];


        %dt1=0.01;
        %timevect=(data1(1,1):dt1:data1(end,1))*DT;
        %timevect2=timevect2(1):dt1:timevect2(end);

        Rss=[];


        figure(10)
        subplot(1,length(topmodels1),cc1)

        % generate forecast curves from each bootstrap realization
        for realization=1:M


            rs_hat=Phatss(realization,1:npatches);
            ps_hat=Phatss(realization,npatches+1:2*npatches);
            as_hat=Phatss(realization,2*npatches+1:3*npatches);
            Ks_hat=Phatss(realization,3*npatches+1:4*npatches);


            IC=zeros(npatches,1);


            if method1>=3
                factor1=Phatss(realization,end-1);
            end

            if method1==5
                d_hat=Phatss(realization,end);
            end



            if onset_thr>0
                IC(1,1)=data1(1,2);
                IC(2:end,1)=1;

                invasions=zeros(npatches,1);
                timeinvasions=zeros(npatches,1);
                Cinvasions=zeros(npatches,1);

                invasions(1)=1;
                timeinvasions(1)=0;
                Cinvasions(1)=0;
            else
                IC(1:end,1)=data1(1,2)./length(IC(1:end,1));

                invasions=zeros(npatches,1);
                timeinvasions=zeros(npatches,1);
                Cinvasions=zeros(npatches,1);

                invasions(1:end)=1;
                timeinvasions(1:end)=0;
                Cinvasions(1:end)=0;
            end



            [~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect2,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);


            for j=1:npatches

                incidence1=[x(1,j);diff(x(:,j))];

                plot(timevect2,incidence1,color1(j,:))
                hold on

            end

            y=sum(x,2);

            totinc=[y(1,1);diff(y(:,1))];

            totinc(1)=totinc(1)-(npatches-1);

            bestfit=totinc;

            gray1=gray(10);


            plot(timevect2,totinc,'color',gray1(7,:))
            %plot(timevect,data1(:,2),'ko')

            curvesforecasts1=[curvesforecasts1 totinc];

            forecasts2=AddPoissonError(cumsum(totinc),20,dist1,factor1,d_hat);

            curvesforecasts2=[curvesforecasts2 forecasts2];


            % <=========================================================================================>
            % <================================ Estimate R_t ===========================================>
            % <=========================================================================================>

            [Rs,ts]=get_Rt(timevect2,totinc,type_GId1,mean_GI1,var_GI1);

            Rss=[Rss Rs];

        end


        % <=============================================================================================>
        % <================================ Plot short-term forecast ====================================>
        % <==============================================================================================>



        line1=plot(data1(:,1),data1(:,2),'ko')
        set(line1,'LineWidth',2)



        axis([0 length(timevect2)-1 0 max(data1(:,2))*2])


        line2=[timevect(end) 0;timevect(end) max(data1(:,2))*2];


        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)


        %caddate1=caddate1(6:end);

        datenum1=datenum([str2num(caddate1(7:10)) str2num(caddate1(1:2)) str2num(caddate1(4:5))]);

        datevec1=datevec(datenum1+forecastingperiod*DT);

        wave=[datevecfirst1 datevec1(1:3)];



        % plot dates in x axis
        'day='
        datenum1=datenum(wave(1:3))+timelags*DT; % start of fall wave (reference date)
        datestr(datenum1)

        datenumIni=datenum1;
        datenumEnd=datenum(wave(4:6))

        dates1=datestr(datenumIni:DT:datenumEnd,'mm-dd');

        if DT==1

            set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));
        elseif DT==7


            set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));

        elseif DT==365

             years1=wave(1)+timelags:wave(4);

            set(gca,'XTick',0:1:length(years1)-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));


        end


        xticklabel_rotate;


        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)

        ylabel(strcat(caddisease,{' '},datatype))

        %title(strcat('Sub-epidemic Model Forecasts',{' '},'-',{' '},datatype,{' '},' - ', {' '},getUSstateName(outbreakx),{' '},'- Reported by ',{' '},caddate1))

        title(strcat(num2ordinal(rank1),' Ranked Model'))

        set(gca,'FontSize',24)
        set(gcf,'color','white')

        % <=========================================================================================>
        % <================================ Plot R_t ===============================================>
        % <=========================================================================================>

        Rtmedian=quantile(Rss',0.5);
        RtLB=quantile(Rss',0.025);
        RtUB=quantile(Rss',0.975);

        figure(502)

        subplot(1,length(topmodels1),cc1)

        plot(ts(2:end),Rss,'c-')
        hold on

        line1=plot(ts(2:end),Rtmedian,'r-')
        set(line1,'LineWidth',2)
        hold on

        line1=plot(ts(2:end),RtLB,'r--')
        set(line1,'LineWidth',2)
        line1=plot(ts(2:end),RtUB,'r--')
        set(line1,'LineWidth',2)

        line2=[0 1;ts(end) 1];
        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)

        axis([0 length(timevect2)-1 0 max(data1(:,2))*2])
        line2=[timevect(end) 0;timevect(end) max(data1(:,2))*2];
        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)

        ylabel('R_t')

        title(strcat(num2ordinal(rankx),' Ranked Model'))


        axis([0 ts(end) 0 2])

        if DT==1

            set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));

        elseif DT==7

            set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));

        elseif DT==365

             years1=wave(1)+timelags:wave(4);

            set(gca,'XTick',0:1:length(years1)-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));

        end


        xticklabel_rotate;

        set(gca,'FontSize', 24);
        set(gcf,'color','white')


        % <=============================================================================================>
        % <============================== Save file with R_t estimates ======================================>
        % <=============================================================================================>

        Rtdata=[ts(2:100:end)' Rtmedian(2:100:end)' RtLB(2:100:end)' RtUB(2:100:end)'];

        T = array2table(Rtdata);
        T.Properties.VariableNames(1:4) = {'time','Rt median','Rt 95%CI LB','Rt 95% CI UB'};
        writetable(T,strcat('./output/Rt-ranked(', num2str(rank1),')-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-area-',num2str(outbreakx),'-',caddate1,'.csv'))


        % <=========================================================================================>
        % <================================ Most recent R estimate =======================>
        % <=========================================================================================>

        param_Rt=[param_Rt;[ts(end) median(Rss(end,:)) quantile(Rss(end,:),0.025) quantile(Rss(end,:),0.975) std(Rss(end,:))]];


        % <=========================================================================================>
        % <================================ Plot the short-term forecast ===========================>
        % <=========================================================================================>

        LB1=quantile(curvesforecasts2',0.025);
        LB1=(LB1>=0).*LB1;

        UB1=quantile(curvesforecasts2',0.975);
        UB1=(UB1>=0).*UB1;


        figure(11)
        subplot(1,length(topmodels1),cc1)

        hold on

        h=area(timevect2',[LB1' UB1'-LB1'],'LineStyle','--')
        hold on

        h(1).FaceColor = [1 1 1];
        h(2).FaceColor = [0.8 0.8 0.8];

        %line1=plot(timevect2,quantile(curvesforecasts2',0.5),'r-')

        line1=plot(timevect2,median(curvesforecasts2,2),'r-')

        set(line1,'LineWidth',2)

        line1=plot(timevect2,LB1,'k--')
        set(line1,'LineWidth',2)

        line1=plot(timevect2,UB1,'k--')
        set(line1,'LineWidth',2)


        gray1=gray(10);

        % plot time series datalatest
        line1=plot(data1(:,1),data1(:,2),'ko')
        set(line1,'LineWidth',2)

        axis([0 length(timevect2)-1 0 max(UB1)*1.2])
        line2=[timevect(end) 0;timevect(end) max(UB1)*1.20];

        box on

        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)

        % plot dates in x axis
        'day='
        datenum1=datenum(wave(1:3))+timelags*DT; % start of fall wave (reference date)
        datestr(datenum1)

        datenumIni=datenum1;
        datenumEnd=datenum(wave(4:6))

        dates1=datestr(datenumIni:DT:datenumEnd,'mm-dd');

        if DT==1

            set(gca, 'XTick', 0:3:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:3:end,:)));

        elseif DT==7

            set(gca, 'XTick', 0:2:length(dates1(:,1))-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',dates1(1:2:end,:)));
            
        elseif DT==365

             years1=wave(1)+timelags:wave(4);

            set(gca,'XTick',0:1:length(years1)-1);
            set(gca, 'XTickLabel', strcat('\fontsize{14}',num2str(years1')));

        end

        xticklabel_rotate;

        line1=plot(line2(:,1),line2(:,2),'k--')
        set(line1,'LineWidth',2)

        ylabel(strcat(caddisease,{' '},datatype))

        %title(strcat('Sub-epidemic Model Forecasts',{' '},'-',{' '},datatype,{' '},' - ', {' '},getUSstateName(outbreakx),{' '},'- Reported by ',{' '},caddate1))

        title(strcat(num2ordinal(rank1),' Ranked Model'))


        set(gca,'FontSize',24)
        set(gcf,'color','white')

        cc1=cc1+1;

    end

end