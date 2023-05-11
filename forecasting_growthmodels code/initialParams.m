
function params0=initialParams(data1,flag1)

switch flag1

    case -1
        params0=[1 1 1 sum(data1) data1(1) 1 1];
        
    case 0
        params0=[1 0.9 1 sum(data1) data1(1) 1 1];

    case 1
        params0=[0.1 0.9 1 sum(data1) data1(1) 1 1];

    case 2
        params0=[0.1 0.9 1 sum(data1) data1(1) 1 1];

    case 3
        params0=[0.1 1 1 sum(data1) data1(1) 1 1];

    case 4
        params0=[0.1 1 0.9 sum(data1) data1(1) 1 1];

    case 5

        %params0=[0.1 1 1 sum(data1) 1];

        K_guess=sum(data1);

        r_guess=1-data1(1)/K_guess;  % r=1-C0/K

        a_guess=r_guess/log(K_guess/data1(1)); % r/log(K/C0)

        params0=[r_guess+0.001 1 a_guess sum(data1) data1(1) 1 1];

    otherwise
        params0=[0.1 0.9 1 sum(data1) data1(1) 1 1];

end

