function curves=AddErrorStructure(yi,M,dist1,factor1,d)

%yi is cumulative curve

curves=[];

for real=1:M
    
    yirData=zeros(length(yi),1);
    
    yirData(1)=yi(1);
    
    switch dist1
      
        case 0
            
            %Normal distribution
            for t=2:length(yi)
                lambda=normrnd(yi(t)-yi(t-1),factor1);
                yirData(t,1)=lambda;
            end
            
        case 1
            
            %Poisson dist
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                yirData(t,1)=poissrnd(lambda,1,1);
            end
            
        case 2
            % Negative binomial dist with VAR=factor*mean1;
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1*factor1;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
        case 3
            % Negative binomial dist with parameter VAR= MEAN + alpha*MEAN
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1+mean1*factor1;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
        case 4
            % Negative binomial dist with parameter VAR= MEAN +
            % alpha*MEAN^2
            
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1+factor1*mean1^2;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
        case 5
            % Negative binomial dist with parameter VAR= MEAN +
            % alpha*MEAN^d
            
            for t=2:length(yi)
                lambda=abs(yi(t)-yi(t-1));
                mean1=lambda;
                var1=mean1+factor1*mean1^d;
                p1=mean1/var1;
                r1=mean1*p1/(1-p1);
                yirData(t,1)=nbinrnd(r1,p1,1,1);
            end
            
    end
    
    curves=[curves yirData];
    
end
