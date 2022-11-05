function [quantilesc1,quantilesf1, alphaquantiles]=computeWIS(data1,curvesforecasts2,forecastingperiod)

%What is Alpha quantile?
%More precisely, an α-quantile is a real number Xα such that the random variable is less or equal to it with probability at least α, and greater or equal to it with probability at least (1–α). It is not possible to require equality because the probability of the value Xα may be positive.May 18, 2011


% calibration period

tf2=length(data1(:,1));

dataCperiod=1:tf2;

alphaquantiles=[0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990];

quantilesc1=zeros(length(dataCperiod),length(alphaquantiles));

for i=1:length(dataCperiod)
    
 
    quantilesc1(i,:)=quantile(curvesforecasts2(i,:),alphaquantiles);
    quantilesc1(i,:)=(quantilesc1(i,:)>=0).*quantilesc1(i,:);

end 


%Forecasting period

quantilesf1=zeros(forecastingperiod,length(alphaquantiles));

if forecastingperiod>0
    
    for forecastingperiod1=forecastingperiod
        
        dataFperiod=tf2+1:tf2+forecastingperiod1;
 
        for i=1:length(dataFperiod)

           quantilesf1(i,:)=quantile(curvesforecasts2(dataFperiod(i),:),alphaquantiles);
           quantilesf1(i,:)=(quantilesf1(i,:)>=0).*quantilesf1(i,:);


        end
        
    end
    
end
