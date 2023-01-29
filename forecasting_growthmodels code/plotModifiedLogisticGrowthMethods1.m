function objfunction=plotModifiedLogisticGrowth1(z)

global flag1 method1 timevect ydata

r=z(1);
p=z(2);
a=z(3);
K=z(4);
I0=z(5);
alpha=z(6);
d=z(7);

%[r p a K I0]

IC=zeros(1,1);

IC(1,1)=I0;

[t,x]=ode15s(@modifiedLogisticGrowth,timevect,IC,[],r,p,a,K,flag1);
incidence1=[x(1,1);diff(x(:,1))];

yfit=incidence1;

eps=0.001;

%%MLE expression
%This is the negative log likelihood, name is legacy from least squares code
%Note that a term that is not a function of the params has been excluded so to get the actual
%negative log-likliehood value you would add: sum(log(factorial(sum(casedata,2))))

if sum(yfit)==0
    objfunction=10^10;%inf;
else
    %    z
    yfit(yfit==0)=eps; %set zeros to eps to allow calculation below.  Shouldn't affect solution, just keep algorithm going.
    
    
    switch method1
        
        case 0  %Least squares
            
            
            %SSE=sum((ydata-yfit).^2);
            
            %objfunction= -1*((-length(ydata)/2)*log(2*pi)-(length(ydata)/2)*log(SSE/length(ydata))-length(ydata)/2);
            
            objfunction=sum((ydata-yfit).^2);
            
            
        case 1 %MLE Poisson (negative log-likelihood)
            
            
            %             sum1=0;
            %             for i=1:length(ydata)
            %
            %                 sum1=sum1+ydata(i)*log(yfit(i))-sum(log(2:1:ydata(i)))-yfit(i);
            %
            %             end
            %
            %             objfunction=-sum1;
            
            objfunction=-sum(ydata.*log(yfit)-yfit);
            
            
        case 3  % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean;
            
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha)*yfit(i));
                    
                end
                
                sum1=sum1+ydata(i)*log(alpha)-(ydata(i)+(1/alpha)*yfit(i))*log(1+alpha)-sum(log(2:1:ydata(i)));
                
            end
            
            objfunction=-sum1;
            
        case 4
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^2;
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha));
                    
                end
                
                sum1=sum1+ydata(i)*log(alpha*yfit(i))-(ydata(i)+(1/alpha))*log(1+alpha*yfit(i))-sum(log(2:1:ydata(i)));
                
            end
            
            objfunction=-sum1;
            
        case 5
            % MLE Negative binomial (negative log-likelihood) where sigma^2=mean+alpha*mean^d;
            
            sum1=0;
            
            for i=1:length(ydata)
                for j=0:(ydata(i)-1)
                    
                    sum1=sum1+log(j+(1/alpha)*yfit(i).^(2-d));
                    
                end
                
                sum1=sum1+ydata(i)*log(alpha*(yfit(i).^(d-1)))-(ydata(i)+(1/alpha)*yfit(i).^(2-d))*log(1+alpha*(yfit(i).^(d-1)))-sum(log(2:1:ydata(i)));
                
            end
            
            objfunction=-sum1;
            
    end
    
    
    %         if ~isreal(objfunction)
    %             dbstop
    %         end
end
