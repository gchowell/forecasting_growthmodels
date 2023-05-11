
function [ratios,timebins1]=getMeanVarianceRatio(data1,binsize1,flag1)

ratios=[];

timebins1=[];

switch flag1
    
    case 1 % no moving window
        
        for i=1:binsize1:length(data1)-binsize1
            
            mean1=mean(data1(i:i+binsize1-1));
            var1=var(data1(i:i+binsize1-1));
            
            ratio1=var1/mean1;
            
            ratios=[ratios;[ratio1 mean1 var1]];
            
            timebins1=[timebins1;(i-1)+(binsize1+1)/2];
            
        end
        
        
    case 2 %moving window
        
        for i=1:1:length(data1)-binsize1
            
            mean1=mean(data1(i:i+binsize1-1));
            var1=var(data1(i:i+binsize1-1));
            
            ratio1=var1/mean1;
            
            ratios=[ratios;[ratio1 mean1 var1]];
            
            timebins1=[timebins1;i];
            
        end
        
        
        
end

