function [CM]=Chaoes_Map(Type,i,xk,xi)
CM(i)=xi;
CM(k)=xk;
switch Type
    case 'CM1' %%CHebyshev
        CM(i+1)=cos(i*acos(CM(i)));
    case 'CM2' %%Circle 
        a=0.5;
        b=0.2;
        CM(i+1)=mod(CM(i) + b - (a/(2*pi)) * sin(2*pi*CM(k)), 1);
    case 'CM3' %%Gauss/Mouse
        if CM(i) == 0
            CM(i+1) = 1;
        else
            CM(i+1) = 1 / mod(CM(i), 1);
        end
    case 'CM4' %%Iterative,a∈（0，1）    
        CM(i+1) = sin(a * pi / CM(i));
    case 'CM5'  %%Logistic,x0∈（0，1）
        CM(i+1) = a * CM(i) * (1 - CM(i));
    case 'CM6'  %%  Piecewise,P∈（0，0.5）,x∈（0，1）,P≠0
        if CM(i) >= 0 && CM(i) < p
            CM(i+1) = CM(i) / p;
        elseif CM(i) >= p && CM(i) < 0.5
            CM(i+1) = (CM(i) - p) / (0.5 - p);
        elseif CM(i) >= 0.5 && CM(i) < (1 - p)
            CM(i+1) = (1 - p - CM(i)) / (0.5 - p);
        elseif CM(i) >= (1 - p) && CM(i) < 1
            CM(i+1) = (1 - CM(i)) / p;
        end
    case 'CM7'  %%  Sine,a∈（0，4）
        CM(i+1) = a/4*sin(pi * CM(i));
    case 'CM8'  %%  Singer,mu∈（0.9，1.08）    
        CM(i+1) = mu * (7.68*CM(i) - 23.31*CM(i)^2 + 28.75*CM(i)^3 - 13.302875*CM(i)^4);
    case 'CM9'  %% Sinusoidal
        a=2.3;
        CM(i+1) = a * CM(i)^2 * sin(pi * CM(i));
    case 'CM10'  %%  Tent
        if CM(i) < 0.7
            CM(i+1) = CM(i) / 0.7;
        else
            CM(i+1) = (10/3) * (1 - CM(i));
        end
end

end