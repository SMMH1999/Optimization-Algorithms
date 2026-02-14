function [bestFitness, bestPosition, convergenceCurve] = LOA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Lion Optimization Algorithm (LOA) - Basic Version
    % =========================================================================
    % Description:
    % A population-based attraction-driven optimizer. Individuals are ranked
    % according to attraction values derived directly from fitness.
    % Top 10% move toward the global best, others move toward random members.
    %
    % INPUTS:
    % lb, ub   : Lower/Upper bounds (scalar or 1×dim vector)
    % dim      : Problem dimension
    % nPop     : Population size
    % maxItr   : Maximum iterations
    % objFun   : Objective function handle (minimization)
    %
    % OUTPUTS:
    % bestFitness      : Best fitness found
    % bestPosition     : Best solution vector
    % convergenceCurve : Best fitness per iteration
    % =========================================================================

    %% Initialization
    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end

    eliteRatio = 0.1;

    pop = rand(nPop, dim) .* (ub - lb) + lb;
    fitness = zeros(nPop,1);

    for i = 1:nPop
        fitness(i) = objFun(pop(i,:));
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = pop(idx,:);

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for t = 1:maxItr

        attraction = 1 ./ (1 + fitness);
        [~, sortedIdx] = sort(attraction, 'descend');

        nElite = max(1, round(eliteRatio * nPop));

        for i = 1:nPop
            id = sortedIdx(i);

            if i <= nElite
                pop(id,:) = pop(id,:) + rand(1,dim) .* (bestPosition - pop(id,:));
            else
                j = randi(nPop);
                pop(id,:) = pop(id,:) + rand(1,dim) .* (pop(j,:) - pop(id,:));
            end

            pop(id,:) = max(pop(id,:), lb);
            pop(id,:) = min(pop(id,:), ub);

            fitness(id) = objFun(pop(id,:));
        end

        [currentBest, idx] = min(fitness);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = pop(idx,:);
        end

        convergenceCurve(t) = bestFitness;
    end

end
