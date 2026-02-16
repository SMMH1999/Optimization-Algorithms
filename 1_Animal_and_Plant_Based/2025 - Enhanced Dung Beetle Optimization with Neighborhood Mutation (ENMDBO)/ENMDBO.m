function [bestFitness, bestPosition, convergenceCurve] = ENMDBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Enhanced Dung Beetle Optimization with Neighborhood Mutation (ENMDBO)
    % =========================================================================
    % Author: Mingyang Yu (2025)
    % Paper:
    % Yu, M., et al. (2025). Improved Coverage and Redundancy Management in WSN
    % Using ENMDBO. IEEE Transactions on Network and Service Management.
    %
    % Benchmark Refactored Version (Single-Function Format)
    %
    % INPUTS:
    % lb        : Lower bound (scalar or vector)
    % ub        : Upper bound (scalar or vector)
    % dim       : Dimension
    % nPop      : Population size
    % maxItr    : Maximum iterations
    % objFun    : Objective function handle
    %
    % OUTPUTS:
    % bestFitness
    % bestPosition
    % convergenceCurve
    % =========================================================================

    %% Parameters
    P_percent = 0.2;
    pNum = round(nPop * P_percent);
    max_no_improvement = 10;
    no_improvement_count = 0;

    if numel(lb) == 1, lb = lb * ones(1, dim); end
    if numel(ub) == 1, ub = ub * ones(1, dim); end

    %% Initialization
    x = lb + (ub - lb) .* rand(nPop, dim);
    fit = zeros(nPop,1);

    for i = 1:nPop
        fit(i) = objFun(x(i,:));
    end

    pX = x;
    pFit = fit;
    XX = pX;

    [bestFitness, idx] = min(fit);
    bestPosition = x(idx,:);

    convergenceCurve = zeros(1,maxItr);

    %% ================= Main Loop =================
    for t = 1:maxItr

        [~, worstIdx] = max(fit);
        worse = x(worstIdx,:);

        % -------- Producers (ECST) --------
        for i = 1:pNum

            indexj = randi(nPop);
            aVec = pX(i,:) - bestPosition;
            bVec = pX(indexj,:) - bestPosition;
            denom = norm(aVec)*norm(bVec);
            if denom == 0
                cosValue = 0;
            else
                cosValue = sum(aVec.*bVec)/denom;
            end

            if cosValue < rand
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

        Xnew1  = max(min(bestXX*(1-R), ub), lb);
        Xnew2  = max(min(bestXX*(1+R), ub), lb);
        Xnew11 = max(min(bestPosition*(1-R), ub), lb);
        Xnew22 = max(min(bestPosition*(1+R), ub), lb);

        % -------- Brood Balls --------
        for i = (pNum+1):min(12,nPop)
            x(i,:) = bestXX + rand(1,dim).*(pX(i,:)-Xnew1) + ...
                rand(1,dim).*(pX(i,:)-Xnew2);
            x(i,:) = max(min(x(i,:), Xnew2), Xnew1);
            fit(i) = objFun(x(i,:));
        end

        % -------- Larvae (NSMS) --------
        for i = 13:min(19,nPop)
            if rand < 0.5
                x(i,:) = pX(i,:) + randn.*(pX(i,:) - Xnew11) + ...
                    rand(1,dim).*(pX(i,:) - Xnew22);
            else
                [~, sortedIdx] = sort(fit);
                pos = find(sortedIdx == i);
                nearPos = max(pos-1,1);
                nearIdx = sortedIdx(nearPos);
                x(i,:) = bestPosition - rand*(bestPosition - pX(nearIdx,:));
            end

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

        % -------- Global Update --------
        [currentBest, idx] = min(fit);

        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = x(idx,:);
            no_improvement_count = 0;
        else
            no_improvement_count = no_improvement_count + 1;
        end

        % -------- TTDM --------
        if no_improvement_count >= max_no_improvement
            for i = 1:nPop
                step = (t/maxItr^rand).*randn(1,dim);
                mutation = step .* x(i,:);
                x(i,:) = bestPosition + rand * mutation;
                x(i,:) = max(min(x(i,:), ub), lb);
                fit(i) = objFun(x(i,:));
            end
        end

        % -------- Personal Best --------
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
