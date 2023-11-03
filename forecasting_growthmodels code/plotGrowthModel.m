function [curve,params]=plotGrowthModel(flag1_pass,r_pass,p_pass,a_pass,K_pass,I0_pass,windowsize1_pass)

% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

% Plot growth model solution using the parameter values supplied by the
% user

% <==============================================================================>
% <============================== Growth model =====================================>
% <==============================================================================>

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

if exist('flag1_pass','var')==1 & isempty(flag1_pass)==0
    flag1=flag1_pass;
else
    'flag1 value is required'
    return
end

% <==================================================================================>
% <========================== Parameter values ======================================>
% <==================================================================================>

a1=1;
p1=1;
K1=1;

switch flag1

    case 0

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required for the GGM model'
            return
        end

        if exist('p_pass','var')==1 & isempty(p_pass)==0
            p1=p_pass;
        else
            'p value is required for the GGM model'
            return

        end

        model_name1='Generalized Growth Model';

        params=[r1 p1];

    case 1

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required for the GLM model'
            return
        end

        if exist('p_pass','var')==1 & isempty(p_pass)==0
            p1=p_pass;
        else
            'p value is required for the GLM model'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required for the GLM model'
            return

        end

        model_name1='Generalized Logistic Growth Model';
        
        params=[r1 p1 K1];

    case 2

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required for the Generalized Richards model'
            return
        end

        if exist('p_pass','var')==1 & isempty(p_pass)==0
            p1=p_pass;
        else
            'p value is required for the GLM model'
            return
        end

        if exist('a_pass','var')==1 & isempty(a_pass)==0
            a1=a_pass;
        else
            'a value is required for the Richards model'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required for the Generalized Richards model'
            return
        end

        model_name1='Generalized Richards Model';

        params=[r1 p1 a1 K1];


    case 3

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required for the LM model'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required for the LM model'
            return

        end

        model_name1='Logistic Growth Model';

         params=[r1 K1];

    case 4

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required for the Richards model'
            return
        end

        if exist('a_pass','var')==1 & isempty(a_pass)==0
            a1=a_pass;
        else
            'a value is required for the Richards model'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required for the Richards model'
            return

        end

      model_name1='Richards Model';

      params=[r1 a1 K1];

    case 5 % Gompertz

        if exist('r_pass','var')==1 & isempty(r_pass)==0
            r1=r_pass;
        else
            'r value is required for the Richards model'
            return
        end

        if exist('a_pass','var')==1 & isempty(a_pass)==0
            a1=a_pass;
        else
            'a value is required for the Richards model'
            return
        end

        if exist('K_pass','var')==1 & isempty(K_pass)==0
            K1=K_pass;
        else
            'K value is required for the Richards model'
            return

        end

      model_name1='Gompertz Model';

      params=[r1 a1 K1];


end

% <==================================================================================>
% <========================== Initial condition C(0) ================================>
% <==================================================================================>

if exist('I0_pass','var')==1 & isempty(I0_pass)==0
    I0=I0_pass;
else
    'Initial value C(0) is required'
    return
end

% <==================================================================================>
% <========================== Time span =============================================>
% <==================================================================================>

if exist('windowsize1_pass','var')==1 & isempty('windowsize1_pass')==0
    windowsize1=windowsize1_pass;
else
    'Simulation duration is required'
    return
end

timevect=0:1:windowsize1;

options=[];

[~,F]=ode15s(@modifiedLogisticGrowth,timevect,I0,options,r1,p1,a1,K1,flag1);

curve=abs([F(1,1);diff(F(:,1))]);

line1=plot(timevect,curve,'b-');
set(line1,'LineWidth',2)

xlabel('Time')

ylabel('dC(t)/dt')

title(model_name1)

set(gca,'FontSize', 24);
set(gcf,'color','white')

