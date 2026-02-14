function [bestFitness, bestPosition, convergenceCurve] = QIO(lb, ub, dim, nPop, maxItr, objFun)
    %% Quadratic Interpolation Optimization (QIO)
    %--------------------------------------------------------------------------
    % Author: Zhao et al., 2023
    % Paper: Quadratic Interpolation Optimization (QIO): A new optimization
    %        algorithm based on generalized quadratic interpolation and its
    %        applications to real-world engineering problems.
    %        Computer Methods in Applied Mechanics and Engineering, 116446
    %--------------------------------------------------------------------------
    % INPUTS:
    %   lb        : lower bound (scalar or 1 x dim vector)
    %   ub        : upper bound (scalar or 1 x dim vector)
    %   dim       : problem dimensionality
    %   nPop      : population size
    %   maxItr    : maximum number of iterations
    %   objFun    : handle to objective function, e.g. @(x) sum(x.^2)
    %
    % OUTPUTS:
    %   bestFitness      : best fitness value found
    %   bestPosition     : position corresponding to bestFitness
    %   convergenceCurve : best fitness at each iteration
    %--------------------------------------------------------------------------
    %% Initialization
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end

    popPos = rand(nPop,dim).*(ub-lb) + lb;
    popFit = zeros(nPop,1);
    for i = 1:nPop
        popFit(i) = objFun(popPos(i,:));
    end

    [bestFitness, idx] = min(popFit);
    bestPosition = popPos(idx,:);
    convergenceCurve = zeros(maxItr,1);

    %% Main Loop
    for it = 1:maxItr
        for i = 1:nPop
            newPos = zeros(1,dim);
            if rand > 0.5
                % --- Exploration ---
                K = [1:i-1 i+1:nPop];
                RandInd = randperm(nPop-1,3);
                K1 = K(RandInd(1));
                K2 = K(RandInd(2));
                K3 = K(RandInd(3));
                f1 = popFit(i); f2 = popFit(K1); f3 = popFit(K2);

                for j = 1:dim
                    newPos(j) = GQI(popPos(i,j), popPos(K1,j), popPos(K2,j), ...
                        f1, f2, f3, lb(j), ub(j));
                end

                a = cos(pi/2*it/maxItr);
                b = 0.7*a + 0.15*a*(cos(5*pi*it/maxItr)+1);
                w1 = 3*b*randn;
                newPos = newPos + w1.*(popPos(K3,:) - newPos) + round(0.5*(0.05+rand))*(log(rand/(rand)));

            else
                % --- Exploitation ---
                K = [1:i-1 i+1:nPop];
                RandInd = randperm(nPop-1,2);
                K1 = K(RandInd(1));
                K2 = K(RandInd(2));
                f1 = popFit(K1); f2 = popFit(K2); f3 = bestFitness;

                for j = 1:dim
                    newPos(j) = GQI(popPos(K1,j), popPos(K2,j), bestPosition(j), ...
                        f1, f2, f3, lb(j), ub(j));
                end

                w2 = 3*(1-(it-1)/maxItr)*randn;
                rD = randi(dim);
                newPos = newPos + w2*(bestPosition - round(1+rand)*(ub-lb)./(ub(rD)-lb(rD))*popPos(i,rD));
            end

            newPos = SpaceBound(newPos, ub, lb);
            newFit = objFun(newPos);

            if newFit < popFit(i)
                popFit(i) = newFit;
                popPos(i,:) = newPos;
            end
        end

        % Update global best
        [currentBestFit, idx] = min(popFit);
        if currentBestFit < bestFitness
            bestFitness = currentBestFit;
            bestPosition = popPos(idx,:);
        end
        convergenceCurve(it) = bestFitness;
    end
end

%% --- Generalized Quadratic Interpolation ---
function L = GQI(a,b,c,fa,fb,fc,low,up)
    fabc = [fa fb fc];
    [fijk, ind] = sort(fabc);
    fi = fijk(1); fj = fijk(2); fk = fijk(3);
    ai = ind(1); bi = ind(2); ci = ind(3);
    x = [a b c];
    xi = x(ai); xj = x(bi); xk = x(ci);

    if (xk>=xi && xi>=xj) || (xj>=xi && xi>=xk)
        L = Interpolation(xi,xj,xk,fi,fj,fk,low,up);
    elseif (xk>=xj && xj>=xi)
        I = Interpolation(xi,xj,xk,fi,fj,fk,low,up);
        if I < xj
            L = I;
        else
            L = Interpolation(xi,xj,3*xi-2*xj,fi,fj,fk,low,up);
        end
    elseif (xi>=xj && xj>=xk)
        I = Interpolation(xi,xj,xk,fi,fj,fk,low,up);
        if I > xj
            L = I;
        else
            L = Interpolation(xi,xj,3*xi-2*xj,fi,fj,fk,low,up);
        end
    elseif (xj>=xk && xk>=xi)
        L = Interpolation(xi,2*xi-xk,xk,fi,fj,fk,low,up);
    elseif (xi>=xk && xk>=xj)
        L = Interpolation(xi,2*xi-xk,xk,fi,fj,fk,low,up);
    end
end

%% --- Quadratic Interpolation ---
function L_xmin = Interpolation(xi,xj,xk,fi,fj,fk,l,u)
    a = (xj^2 - xk^2)*fi + (xk^2 - xi^2)*fj + (xi^2 - xj^2)*fk;
    b = 2*((xj - xk)*fi + (xk - xi)*fj + (xi - xj)*fk);
    L_xmin = a/(b+eps);
    if isnan(L_xmin) || isinf(L_xmin) || L_xmin>u || L_xmin<l
        L_xmin = rand*(u-l) + l;
    end
end

%% --- Boundary Control ---
function X = SpaceBound(X,Up,Low)
    Dim = length(X);
    S = (X>Up) + (X<Low);
    X = (rand(1,Dim).*(Up-Low)+Low).*S + X.*(~S);
end
