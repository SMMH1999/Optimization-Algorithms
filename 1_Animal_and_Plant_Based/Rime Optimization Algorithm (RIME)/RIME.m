function [bestFitness, bestPosition, convergenceCurve] = RIME(lb, ub, dim, nPop, maxItr, objFun)
    %====================================================================
    % RIME: Randomized Inspired Metaheuristic for Engineering optimization
    % Physics-based optimization algorithm
    %
    % Authors: Ali Asghar Heidari, Huiling Chen, Hang Su, Dong Zhao, et al.
    %
    % Description:
    %   This function performs optimization using the RIME algorithm.
    %   Minimization problems are assumed.
    %
    % Inputs:
    %   lb      - (scalar or vector) Lower bound of variables
    %   ub      - (scalar or vector) Upper bound of variables
    %   dim     - (int) Dimensionality of the problem
    %   nPop    - (int) Number of search agents (population size)
    %   maxItr  - (int) Maximum number of iterations
    %   objFun  - (function handle) Objective function to minimize
    %
    % Outputs:
    %   bestFitness      - (double) Best objective function value found
    %   bestPosition     - (1 x dim) Position of the best solution
    %   convergenceCurve - (1 x maxItr) Best fitness value at each iteration
    %
    % Tunable parameters:
    %   W - Soft-Rime control parameter (default 5)
    %====================================================================

    %% Initialization
    bestPosition = zeros(1, dim);
    bestFitness = inf; % Change to -inf for maximization
    W = 5; % Soft-Rime parameter

    % Initialize population
    population = initializePopulation(nPop, dim, lb, ub);

    % Fitness initialization
    fitness = zeros(1, nPop);
    for i = 1:nPop
        fitness(i) = objFun(population(i, :));
        if fitness(i) < bestFitness
            bestFitness = fitness(i);
            bestPosition = population(i, :);
        end
    end

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for t = 1:maxItr
        RimeFactor = (rand-0.5)*2*cos(pi*t/(maxItr/10))*(1 - round(t*W/maxItr)/W);
        E = sqrt(t/maxItr);

        newPopulation = population;
        normalizedFitness = normr(fitness);

        for i = 1:nPop
            for j = 1:dim
                % Soft-Rime search
                if rand < E
                    newPopulation(i,j) = bestPosition(j) + RimeFactor*((ub(j)-lb(j))*rand + lb(j));
                end
                % Hard-Rime puncture
                if rand < normalizedFitness(i)
                    newPopulation(i,j) = bestPosition(j);
                end
            end
        end

        % Boundary handling and greedy selection
        for i = 1:nPop
            newPopulation(i,:) = max(min(newPopulation(i,:), ub), lb);

            newFitness = objFun(newPopulation(i,:));

            if newFitness < fitness(i)
                fitness(i) = newFitness;
                population(i,:) = newPopulation(i,:);
                if newFitness < bestFitness
                    bestFitness = newFitness;
                    bestPosition = newPopulation(i,:);
                end
            end
        end

        convergenceCurve(t) = bestFitness;
    end

end

%% Helper Function: Initialize Population
function positions = initializePopulation(nPop, dim, lb, ub)
    % Supports scalar or vector bounds
    if isscalar(lb)
        positions = rand(nPop, dim).*(ub-lb) + lb;
    else
        positions = zeros(nPop, dim);
        for i = 1:dim
            positions(:, i) = rand(nPop,1).*(ub(i)-lb(i)) + lb(i);
        end
    end
end
