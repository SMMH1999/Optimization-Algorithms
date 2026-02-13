function [bestFitness, bestPosition, convergenceCurve] = DBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Dung Beetle Optimizer (DBO)
    % =========================================================================
    % Original Authors: Jiankai Xue & Bo Shen (2022)
    % Source Paper:
    % Xue, J., & Shen, B. (2022). Dung beetle optimizer: a new meta-heuristic
    % algorithm for global optimization. Journal of Supercomputing.
    %
    % Benchmark Refactored Version (Single-Function Format)
    %
    % INPUTS:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Problem dimension
    % nPop      : Population size
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % OUTPUTS:
    % bestFitness      : Best objective value found
    % bestPosition     : Best solution vector
    % convergenceCurve : Best fitness value at each iteration
    %
    % =========================================================================

    %% -------------------- Parameters --------------------
    P_percent = 0.2;
    pNum = round(nPop * P_percent);

    if numel(lb) == 1, lb = lb * ones(1, dim); end
    if numel(ub) == 1, ub = ub * ones(1, dim); end

    %% -------------------- Initialization --------------------
    x = zeros(nPop, dim);
    fit = zeros(nPop, 1);

    for i = 1:nPop
        x(i,:) = lb + (ub - lb) .* rand(1, dim);
        fit(i) = objFun(x(i,:));
    end

    pX = x;
    pFit = fit;
    XX = pX;

    [bestFitness, idx] = min(fit);
    bestPosition = x(idx,:);

    convergenceCurve = zeros(1, maxItr);

    %% ==================== Main Loop ====================
    for t = 1:maxItr

        [~, worstIdx] = max(fit);
        worse = x(worstIdx,:);
        r2 = rand;

        % -------- Producers --------
        for i = 1:pNum

            if r2 < 0.9
                a = 1;
                if rand <= 0.1
                    a = -1;
                end
                x(i,:) = pX(i,:) + 0.3*abs(pX(i,:) - worse) + a*0.1*XX(i,:);
            else
                aaa = randperm(180,1);
                theta = aaa*pi/180;
                x(i,:) = pX(i,:) + tan(theta).*abs(pX(i,:) - XX(i,:));
            end

            x(i,:) = max(min(x(i,:), ub), lb);
            fit(i) = objFun(x(i,:));
        end

        [~, bestII] = min(fit);
        bestXX = x(bestII,:);

        R = 1 - t/maxItr;

        Xnew1 = max(min(bestXX*(1-R), ub), lb);
        Xnew2 = max(min(bestXX*(1+R), ub), lb);

        Xnew11 = max(min(bestPosition*(1-R), ub), lb);
        Xnew22 = max(min(bestPosition*(1+R), ub), lb);

        % -------- Brood Balls --------
        for i = (pNum+1):min(12,nPop)
            x(i,:) = bestXX + rand(1,dim).*(pX(i,:)-Xnew1) + ...
                rand(1,dim).*(pX(i,:)-Xnew2);
            x(i,:) = max(min(x(i,:), Xnew2), Xnew1);
            fit(i) = objFun(x(i,:));
        end

        % -------- Larvae --------
        for i = 13:min(19,nPop)
            x(i,:) = pX(i,:) + randn.*(pX(i,:) - Xnew11) + ...
                rand(1,dim).*(pX(i,:) - Xnew22);
            x(i,:) = max(min(x(i,:), ub), lb);
            fit(i) = objFun(x(i,:));
        end

        % -------- Thieves --------
        for j = 20:nPop
            x(j,:) = bestPosition + randn(1,dim).* ...
                (abs(pX(j,:) - bestXX) + abs(pX(j,:) - bestPosition))/2;
            x(j,:) = max(min(x(j,:), ub), lb);
            fit(j) = objFun(x(j,:));
        end

        % -------- Update Personal & Global Best --------
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

end
