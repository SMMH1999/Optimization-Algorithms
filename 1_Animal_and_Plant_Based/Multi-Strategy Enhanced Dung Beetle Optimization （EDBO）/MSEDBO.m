function [bestFitness, bestPosition, convergenceCurve] = MSEDBO(lb, ub, dim, nPop, maxItr, objFun)
    % -----------------------------------------------------------------------------------------------------------
    % Multi-Strategy Enhanced Dung Beetle Optimization (EDBO)
    % -----------------------------------------------------------------------------------------------------------
    % Author: Original algorithm author
    %
    % Description:
    %   EDBO is a population-based metaheuristic optimization algorithm inspired
    %   by dung beetle behavior with enhanced exploration and exploitation strategies.
    %
    % Inputs:
    %   lb        - scalar or 1xD vector of lower bounds
    %   ub        - scalar or 1xD vector of upper bounds
    %   dim       - problem dimension
    %   nPop      - population size
    %   maxItr    - maximum number of iterations
    %   objFun    - handle to objective function
    %
    % Outputs:
    %   bestFitness       - best fitness value found
    %   bestPosition      - corresponding decision variable vector
    %   convergenceCurve  - best fitness at each iteration
    %
    % Tunable Parameters:
    %   P_percent - fraction of population considered as "producers" (default 0.2)
    % -----------------------------------------------------------------------------------------------------------

    %% Ensure bounds are vectors
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end

    %% Parameters
    P_percent = 0.2;                  % Fraction of population as producers
    pNum = round(nPop * P_percent);   % Number of producer individuals

    %% Initialization
    x = lb + (ub - lb) .* rand(nPop, dim);  % Random population
    fit = zeros(nPop,1);

    for i = 1:nPop
        fit(i) = objFun(x(i,:));
    end

    pFit = fit;
    pX = x;
    XX = pX;

    [bestFitness, bestI] = min(fit);
    bestPosition = x(bestI,:);

    convergenceCurve = zeros(1,maxItr);

    %% Main loop
    for t = 1:maxItr

        [fmax, B] = max(fit);
        worse = x(B,:);
        r2 = rand;

        %% Producer update
        for i = 1:pNum
            if r2 < 0.9
                a = 1; if rand <= 0.1, a = -1; end
                x(i,:) = pX(i,:) + 0.3*abs(pX(i,:) + rand*bestPosition) + a*0.1*(XX(i,:));
            else
                thetaDeg = randperm(180,1);
                theta = thetaDeg*pi/180;
                x(i,:) = pX(i,:) + tan(theta).*abs(pX(i,:) - XX(i,:));
            end

            % Boundary handling
            x(i,:) = BoundCheck(x(i,:), lb, ub, bestPosition);

            % Fitness evaluation
            fit(i) = objFun(x(i,:));
        end

        [fMMin, bestII] = min(fit);
        bestXX = x(bestII,:);

        R = 1 - rand*(t/maxItr)^2;

        %% Auxiliary positions
        Xnew1 = BoundCheck(bestXX*(1-R), lb, ub);
        Xnew2 = BoundCheck(bestXX*(1+R), lb, ub);
        Xnew11 = BoundCheck(bestPosition*(1-R), lb, ub);
        Xnew22 = BoundCheck(bestPosition*(1+R), lb, ub);

        %% Middle population update
        for i = (pNum+1):12
            x(i,:) = bestXX + rand(1,dim).*(pX(i,:) - Xnew1) + rand(1,dim).*(pX(i,:) - Xnew2);
            x(i,:) = BoundCheck(x(i,:), lb, ub, bestPosition);
            fit(i) = objFun(x(i,:));
        end

        %% Advanced update
        for i = 13:19
            if rand < 0.5
                x(i,:) = pX(i,:) + randn*(pX(i,:) - Xnew11) + rand(1,dim).*(pX(i,:) - Xnew22);
            else
                I = round(1 + rand);
                x(i,:) = pX(i,:) + rand*(bestXX - I*pX(i,:));
            end
            x(i,:) = BoundCheck(x(i,:), lb, ub, bestPosition);
            fit(i) = objFun(x(i,:));
        end

        %% Remaining population update
        for i = 20:nPop
            x(i,:) = bestPosition + randn*(abs(pX(i,:) - bestXX) + abs(pX(i,:) - bestPosition))/2;
            x(i,:) = BoundCheck(x(i,:), lb, ub, bestPosition);
            fit(i) = objFun(x(i,:));
        end

        %% Update personal and global bests
        XX = pX;
        for i = 1:nPop
            if fit(i) < pFit(i)
                pFit(i) = fit(i);
                pX(i,:) = x(i,:);
            end
            if pFit(i) < bestFitness
                bestFitness = pFit(i);
                bestPosition = pX(i,:);
            end
        end

        convergenceCurve(t) = bestFitness;
    end

    %% --- Nested helper function for bounds handling ---
    function s = BoundCheck(s,Lb,Ub,bestPos)
        I = s < Lb; s(I) = rand*bestPos(I);
        J = s > Ub; s(J) = (1-rand)*bestPos(J);
    end

end
