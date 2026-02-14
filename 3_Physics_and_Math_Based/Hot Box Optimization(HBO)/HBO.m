function [bestFitness, bestPosition, convergenceCurve] = Algorithm(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Hot Box Optimization (HBO)
    % Population-based Simulated Annealing Variant
    %
    % Author: Unknown (Refactored for Benchmark Framework)
    %
    % Description:
    % A population-based stochastic optimization algorithm inspired by
    % Simulated Annealing. Each individual is perturbed using a temperature-
    % controlled random displacement and accepted based on Metropolis criterion.
    %
    % =========================================================================
    % INPUTS:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Number of decision variables
    % nPop      : Population size
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % OUTPUTS:
    % bestFitness      : Best fitness value found
    % bestPosition     : Best solution vector found
    % convergenceCurve : Best fitness at each iteration
    %
    % =========================================================================
    % Tunable Parameters:
    % initialTemp  : Initial temperature (default = 100)
    % coolingRate  : Temperature reduction factor (default = 0.99)
    % minTemp      : Minimum stopping temperature (default = 1e-3)
    % =========================================================================

    %% -------------------- Parameter Handling --------------------

    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end

    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    initialTemp = 100;
    coolingRate = 0.99;
    minTemp     = 1e-3;

    %% -------------------- Initialization --------------------

    pop = lb + (ub - lb) .* rand(nPop, dim);
    fitness = zeros(nPop,1);

    for i = 1:nPop
        fitness(i) = objFun(pop(i,:));
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = pop(idx,:);

    convergenceCurve = zeros(maxItr,1);

    %% -------------------- Main Loop --------------------

    for t = 1:maxItr

        temp = initialTemp * (coolingRate^t);

        if temp < minTemp
            convergenceCurve(t:end) = bestFitness;
            break;
        end

        for i = 1:nPop

            % Generate new candidate
            newSol = pop(i,:) + (rand(1,dim) - 0.5) .* temp;

            % Boundary control
            newSol = max(newSol, lb);
            newSol = min(newSol, ub);

            % Evaluate
            newFitness = objFun(newSol);

            % Metropolis acceptance rule
            if newFitness < fitness(i) || rand < exp(-(newFitness - fitness(i))/temp)
                pop(i,:) = newSol;
                fitness(i) = newFitness;
            end
        end

        %% --------- Global Best Update ---------
        [currBest, idx] = min(fitness);

        if currBest < bestFitness
            bestFitness = currBest;
            bestPosition = pop(idx,:);
        end

        %% --------- Convergence Curve ---------
        convergenceCurve(t) = bestFitness;

    end

end
