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

w  = min(7, size(data1,1));                    % early window

if w < 2                                       % degenerate: too few points
    r0  = params0(1);
    rlb = 1e-3;  rub = 10;
else
    y  = max(data1(1:w,2), 0.5);               % guard zeros
    b  = polyfit(data1(1:w,1)*DT, log(y), 1);

    if b(1) <= 0
        r0  = params0(1);      % no growth signal: keep supplied guess
        rlb = 1e-3;  rub = 10; % wide default
    else
        r0  = b(1);
        rlb = r0/20; rub = r0*20;
    end
end

%rlb=median(abs(data1(1:2,2)))/200;
%rub=5*max(abs(data1(1:2,2)));

Cobs = sum(abs(data1(:,2)));
Klb  = max(20, Cobs);            % can't undershoot observed cumulative
Kub  = min(1e3*Cobs, 1e9);       % capped


%Kmax=1000000000;

if fixI0==1

    switch flag1

         case -1  %EXP
            LB=[rlb  1 1 0 I0 LBe];
            UB=[rub  1 1 0 I0 UBe];

        case 0   %GGM
            LB=[rlb  0.01 1 0 I0 LBe];
            UB=[rub  1 1 0 I0 UBe];

        case 1 % GLM
            LB=[rlb  0.01 1 Klb I0 LBe];
            UB=[rub  1 1 Kub I0 UBe];

        case 2 %GRM
            LB=[rlb  0.01 0 Klb I0 LBe];
            UB=[rub  1 10 Kub I0 UBe];

        case 3 %Logistic
            LB=[rlb  1 1 Klb I0 LBe];
            UB=[rub  1 1 Kub I0 UBe];

        case 4 % Richards
            LB=[rlb  1 0 Klb I0 LBe];
            UB=[rub  1 10 Kub I0 UBe];

        case 5 % Gompertz
            LB=[max(1e-4, params0(1)/10)  1 max(1e-6, params0(3)/10) 1 I0 LBe];
            UB=[params0(1)*10             1 params0(3)*10            1 I0 UBe];

    end

else

    y1 = data1(1,2);

    I0lb=max(1, 0.2*y1);
    I0ub=5*max(1,y1);
    I0ub=min(I0ub, Klb);         % keep I0 <= K for saturating models
    I0lb=min(I0lb, I0ub);        % preserve I0lb <= I0ub

    switch flag1

        case -1
            LB=[rlb 1 1 0 I0lb LBe];
            UB=[rub  1 1 0 I0ub UBe];

        case 0
            LB=[rlb 0.01 1 0 I0lb LBe];
            UB=[rub  1 1 0 I0ub UBe];

        case 1
            LB=[rlb  0.01 1 Klb I0lb LBe];
            UB=[rub 1 1 Kub I0ub UBe];

        case 2
            LB=[rlb  0.01 0 Klb I0lb LBe];
            UB=[rub  1 10 Kub I0ub UBe];

        case 3
            LB=[rlb  1 1 Klb I0lb LBe];
            UB=[rub  1 10 Kub I0ub UBe];

        case 4
            LB=[rlb  1 0 Klb I0lb LBe];
            UB=[rub  1 10 Kub I0ub UBe];

        case 5 % Gompertz
            LB=[max(1e-4, params0(1)/10)  1 max(1e-6, params0(3)/10) 1 I0lb LBe];
            UB=[params0(1)*10             1 params0(3)*10            1 I0ub UBe];

    end

end

if flag1 ~= 5
    z(1) = min(max(r0, LB(1)), UB(1));   % informed r guess
end
z = min(max(z, LB), UB);                 % clamp everything into bounds


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

maxAttempts=3; % Initial attempt plus at most two retries.
attempt=0;
flagg=-1;

while flagg<0 && attempt<maxAttempts

    attempt=attempt+1;

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

if flagg<0
    error('GrowthPredict:fit_model:MaxAttemptsReached', ...
        'Fit failed after %d attempts (last exit flag: %g).', ...
        attempt,flagg);
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


