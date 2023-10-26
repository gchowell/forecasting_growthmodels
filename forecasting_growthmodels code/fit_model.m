% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function [P residual fitcurve forecastcurve timevect2,initialguess,fval]=fit_model(data1,params0,fixI0,numstartpoints,DT,flagX,forecastingperiod)

global flag1 method1 timevect ydata

flag1=flagX;

timevect=data1(:,1)*DT;

%timevect=(data1(:,1))*DT;

r=params0(1);
p=params0(2);
a=params0(3);
K=params0(4);
I0=params0(5);
alpha=params0(6);
d=params0(7);

I0=data1(1,2); % initial condition

z(1)=r;
z(2)=p;
z(3)=a;
z(4)=K;
z(5)=I0;
z(6)=alpha;
z(7)=d;


switch method1
    case 0
        LBe=[0 0];
        UBe=[0 0];
    case 1
        LBe=[0 0];
        UBe=[0 0];
    case 3
        LBe=[10^-8 1];
        UBe=[10^3 1];
    case 4
        LBe=[10^-8 1];
        UBe=[10^5 1];
    case 5
        LBe=[10^-8 0.6]; %d>=1
        UBe=[10^5 10^3];

end


rlb=mean(abs(data1(1:2,2)))/200;
rub=mean(abs(data1(1:2,2)))*2;

Kmax=100000000000;

if fixI0==1

    switch flag1

         case -1  %EXP
            LB=[rlb  1 1 0 I0 LBe];
            UB=[rub  1 1 0 I0 UBe];

        case 0   %GGM
            LB=[rlb  0.01 1 0 I0 LBe];
            UB=[rub  1 1 0 I0 UBe];

        case 1 % GLM
            LB=[rlb  0.01 1 20 I0 LBe];
            UB=[rub  1 1 Kmax I0 UBe];

        case 2 %GRM
            LB=[rlb  0.01 0 20 I0 LBe];
            UB=[rub  1 10 Kmax I0 UBe];

        case 3 %Logistic
            LB=[rlb  1 1 20 I0 LBe];
            UB=[rub  1 1 Kmax I0 UBe];

        case 4 % Richards
            LB=[rlb  1 0 20 I0 LBe];
            UB=[rub  1 10 Kmax I0 UBe];

        case 5 % Gompertz
            LB=[rlb  1 0 1 I0 LBe];
            UB=[rub  1 params0(3)+5 1 I0 UBe];

    end

else
    
    switch flag1

        case -1
            LB=[rlb 1 1 0 1 LBe];
            UB=[rub  1 1 0 sum(abs(data1(:,2))) UBe];

        case 0
            LB=[rlb 0.01 1 0 1 LBe];
            UB=[rub  1 1 0 sum(abs(data1(:,2))) UBe];

        case 1
            LB=[rlb  0.01 1 20 1 LBe];
            UB=[rub 1 1 Kmax sum(abs(data1(:,2))) UBe];

        case 2
            LB=[rlb  0.01 0 20 1 LBe];
            UB=[rub  1 10 Kmax sum(abs(data1(:,2))) UBe];

        case 3
            LB=[rlb  1 1 20 1 LBe];
            UB=[rub  1 10 Kmax sum(abs(data1(:,2))) UBe];

        case 4
            LB=[rlb  1 0 20 1 LBe];
            UB=[rub  1 10 Kmax sum(abs(data1(:,2))) UBe];

        case 5
            LB=[rlb  1 0 1 1 LBe];
            UB=[rub  1 params0(3)+5 1 sum(abs(data1(:,2))) UBe];

    end

end


% if 0 % USE LSQCURVEFIT (Non-linear least squares)
% 
%     options=optimset('tolfun',10^-5,'TolX',10^-5,'MaxFunEvals',3200,'MaxIter',3200, 'algorithm','trust-region-reflective');
% 
%     [P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowth1,z,timevect,data1(:,2),LB,UB,options);
% 
%     f=@plotModifiedLogisticGrowth1;
% 
%     problem = createOptimProblem('lsqcurvefit','x0',z,'objective',f,'lb',LB,'ub',UB,'xdata',timevect,'ydata',data1(:,2),'options',options);
% 
%     %ms = MultiStart('PlotFcns',@gsplotbestf,'Display','final');
% 
%     ms = MultiStart('Display','final');
% 
%     ms = MultiStart(ms,'StartPointsToRun','bounds')
% 
%     [P,errormulti] = run(ms,problem,20)
% 
%     z=P;
% 
%     [P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowth1,z,timevect,data1(:,2),LB,UB,options);
% 
% end

%A=[];      % We are using fmincon, but using none of the constraint options
%b=[];
%Aeq=[];
%beq=[];
%nonlcon=[];

%options=optimset('tolfun',10^-5,'TolX',10^-5,'MaxFunEvals',3200,'MaxIter',3200, 'algorithm','interior-point');

%[P, fval, exitflag]=fmincon(@plotModifiedLogisticGrowthMethods1,z,A,b,Aeq,beq,LB,UB,nonlcon,options);

%method1=1; %LSQ=0, MLE (Poisson)=1, Pearson chi-squared=2. MLE(neg binomial)=3

ydata=data1(:,2);

options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-6,'MaxFunEvals',20000,'MaxIter',20000);

%options=optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',10000,'MaxIter',10000);

f=@plotModifiedLogisticGrowthMethods1;

problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);

%ms = MultiStart('PlotFcns',@gsplotbestf);
%ms = MultiStart('Display','final');
ms = MultiStart('Display','off');

%pts = z;
tpoints = CustomStartPointSet(z);

flagg=-1;

while flagg<0

    initialguess=[];

    if numstartpoints>0
        rpoints = RandomStartPointSet('NumStartPoints',numstartpoints); % start with a few random starting sets in addition to the guess supplied by the user (z)

        allpts = {rpoints,tpoints};
        initialguess=list(rpoints,problem);

    else
        allpts = {tpoints};

    end

    initialguess=[initialguess;z];

    %z
    %list(tpoints)

    %ms = MultiStart(ms,'StartPointsToRun','bounds')
    %[xmin,fmin,flag,outpt,allmins] = run(ms,problem,allpts);

    [P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);

end

% ydata
% initialguess
% P
% pause

% P is the vector with the estimated parameters
r_hat=P(1);
p_hat=P(2);
a_hat=P(3);
K_hat=P(4);
I0_hat=P(5);
alpha_hat=P(6);
d_hat=P(7);

options = [];

[~,F]=ode15s(@modifiedLogisticGrowth,timevect,I0_hat,options,r_hat,p_hat,a_hat,K_hat,flag1);

fitcurve=abs([F(1,1);diff(F(:,1))]);

residual=fitcurve-ydata;

%fitcurve=residual+data1(:,2);

if forecastingperiod<1

    forecastcurve=residual+data1(:,2);
    timevect2=timevect;

else

    timevect2=(data1(1,1):data1(end,1)+forecastingperiod)*DT;

    [~,F]=ode15s(@modifiedLogisticGrowth,timevect2,I0_hat,[],r_hat,p_hat,a_hat,K_hat,flag1);

    forecastcurve=abs([F(1,1);diff(F(:,1))]);

end


