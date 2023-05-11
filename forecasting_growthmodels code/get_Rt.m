function [Rs,ts]=get_Rt(timevect,incidence1,type_GId1,mean1,var1)

% type of GI distribution 1=Gamma, 2=Exponential, 3=Delta


dt1=0.01;
timevect2=timevect(1):dt1:timevect(end);

incidence1=interp1(timevect,incidence1,timevect2);

timevect=timevect2;

ts=timevect;


if type_GId1==1
    
    
    b=var1/mean1;
    a=mean1/b;
    
    %a=4.2
    %b=0.68
    prob1=gamcdf(timevect,a,b);
    prob1=diff(prob1)';
    
    sum(prob1);
end

%exponential

if type_GId1==2
    
    
    prob1=expcdf(timevect,mean1);
    
    prob1=diff(prob1)';
    sum(prob1);
    
end


%delta
if type_GId1==3
    
    %mean1=3;
    prob1=zeros(length(timevect),1);
    prob1(mean1/dt1+1)=1;
    sum(prob1);
    
end


Rs=[];

for i=2:length(incidence1)
    
    sum1=0;
    
    for j=1:i-1
        %sum1=sum1+(incidence1(i)-incidence1(j))*prob1(j,1);
        sum1=sum1+(incidence1(i-j))*prob1(j,1);
        %'prob1='
        %prob1(j,1)
        
    end
    
    if sum1>0
        R1=incidence1(i)/sum1;
        
    else
        R1=NaN;
        
    end
    
    Rs=[Rs;R1];
    
end
