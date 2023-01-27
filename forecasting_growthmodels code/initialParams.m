
function params0=initialParams(f_ensembledata,flag1,fixI0)

switch flag1

    case 0
        params0=[1 0.9 1 sum(f_ensembledata) 1 1 1];

    case 1
        params0=[0.1 0.9 1 sum(f_ensembledata) 1 1 1];

    case 2
        params0=[0.1 0.9 1 sum(f_ensembledata) 1 1 1];

    case 3
        params0=[0.1 1 1 sum(f_ensembledata) 1 1 1];

    case 4
        params0=[0.1 1 0.9 sum(f_ensembledata) 1 1 1];

    case 5

        %params0=[0.1 1 1 sum(f_ensembledata) 1];

        K_guess=sum(f_ensembledata);

        r_guess=1-f_ensembledata(1)/K_guess;  % r=1-C0/K

        a_guess=r_guess/log(K_guess/f_ensembledata(1)); % r/log(K/C0)

        params0=[r_guess 1 a_guess sum(f_ensembledata) 1 1 1];

    otherwise
        params0=[0.1 0.9 1 sum(f_ensembledata) 1 1 1];

end

