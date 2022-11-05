function dx=modifiedLogisticGrowth(t,x,r,p,a,K,flag1)


dx=zeros(1,1);

if flag1==0 % GGM
    
    dx(1,1)=r*x(1,1).^p;
    
elseif flag1==1 % GLM
    
    dx(1,1)=r*(x(1,1).^p).*(1-(x(1,1)/K));
    %dx(2,1)=(x(1,1)>K*(1-0.99)).*(r*(x(2,1).^p).*(1-(x(2,1)/K)));
    
    
elseif flag1==2 %GRM
    
    dx(1,1)=r*(x(1,1).^p)*(1-(x(1,1)/K).^a);
    
    
elseif flag1==3 % Logistic
    
    dx(1,1)=r*x.*(1-(x/K));
    
elseif flag1==4 %Richards
    
    dx(1,1)=r*x*(1-(x/K).^a);

elseif flag1==5 % Gompertz
    
    dx(1,1)=r*x*exp(-a*t);
    
    %dx(1,1)=r*x*(log(K/x)).^0.6;
    
elseif flag1==6
    
    
end



% if flag1==0 % GGM
%     
%     dx(1,1)=r*(x.^p);
%     
% elseif flag1==1 % GLM
%     
%     dx(1,1)=r*K*((x/K).^p).*(1-(x/K));
%     
% 
% elseif flag1==2 %GRM
%     
%     dx(1,1)=r*(x.^p)*(1-(x/K).^a);
%     
%     
% elseif flag1==3 % Logistic
%     
%     dx(1,1)=r*x.*(1-(x/K));
%     
% elseif flag1==4 %Richards
%     
%     dx(1,1)=r*x*(1-(x/K).^a);
%     
%     
% end


