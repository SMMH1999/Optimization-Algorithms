function [bestFitness, bestPosition, convergenceCurve] = HSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Holistic Swarm Optimization (HSO)
    % =========================================================================
    % Based on:
    % "Holistic Swarm Optimization: A Novel Metaphor-less Algorithm Guided by
    % Whole Population Information for Addressing Exploration-Exploitation Dilemma"
    % DOI: https://doi.org/10.1016/j.cma.2025.118208
    %
    % -------------------------------------------------------------------------
    % Description:
    % Holistic Swarm Optimization (HSO) is a population-based metaheuristic
    % algorithm that updates each agent using whole-population information.
    % It integrates simulated annealing-based acceptance and adaptive mutation
    % to balance exploration and exploitation.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best decision variable vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   alpha                 - Scaling factor for position updates (default=3)
    %   initialTemp           - Initial temperature for simulated annealing
    %   coolingRate           - Cooling rate for temperature decay
    %   initialMutationRate   - Initial probability of mutation
    %   finalMutationRate     - Final mutation probability
    %   initialMutationStep   - Initial mutation step size
    %   finalMutationStep     - Final mutation step size
    %
    % =========================================================================

    %% Parameter Settings
    alpha = 3;

    % Simulated Annealing parameters
    initialTemp = 10000;
    coolingRate = 0.995;

    % Adaptive Mutation parameters
    initialMutationRate = 0.5;
    finalMutationRate   = 0.1;
    initialMutationStep = 0.3;
    finalMutationStep   = 0.1;

    %% Bound Handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% Initialization
    positions = rand(nPop, dim) .* (ub - lb) + lb;
    positionsNew = positions;

    fitness = inf(1, nPop);
    fitnessNew = inf(1, nPop);

    bestPosition = zeros(1, dim);
    bestFitness = inf;

    convergenceCurve = zeros(1, maxItr);

    % Initial Fitness Evaluation
    for i = 1:nPop
        fitness(i) = objFun(positions(i, :));
    end

    [fitness, idx] = sort(fitness);
    positions = positions(idx, :);

    bestFitness = fitness(1);
    bestPosition = positions(1, :);

    %% Main Loop
    for t = 1:maxItr

        temp = initialTemp * (coolingRate ^ (t-1));

        % Adaptive mutation parameters
        mutationRate = initialMutationRate - (t-1) * ...
            ((initialMutationRate - finalMutationRate) / maxItr);

        mutationStep = initialMutationStep - (t-1) * ...
            ((initialMutationStep - finalMutationStep) / maxItr);

        %% Position Update (Holistic Interaction)
        fitnessVector = fitness(:);
        differences = rms(fitnessVector) - fitnessVector;
        updateCoef = differences / sum(abs(differences) + eps);

        for i = 1:nPop
            for j = 1:dim
                randWeights = rand(nPop, 1);
                displacement = alpha * ...
                    sum(randWeights .* updateCoef .* ...
                    (positions(:, j) - positions(i, j)));

                positionsNew(i, j) = positions(i, j) + displacement;
            end
        end

        %% Adaptive Mutation
        for i = 1:nPop
            if rand < mutationRate
                mutationVector = mutationStep * randn(1, dim);
                positionsNew(i, :) = positionsNew(i, :) + mutationVector;
            end
        end

        %% Boundary Control
        positionsNew = max(positionsNew, lb);
        positionsNew = min(positionsNew, ub);

        %% Fitness Evaluation and SA-based Selection
        for i = 1:nPop
            fitnessNew(i) = objFun(positionsNew(i, :));

            if fitnessNew(i) < fitness(i)
                positions(i, :) = positionsNew(i, :);
                fitness(i) = fitnessNew(i);
            else
                delta = fitnessNew(i) - fitness(i);
                if exp(-delta / temp) > rand
                    positions(i, :) = positionsNew(i, :);
                    fitness(i) = fitnessNew(i);
                end
            end
        end

        %% Sort Population
        [fitness, idx] = sort(fitness);
        positions = positions(idx, :);

        %% Update Global Best
        if fitness(1) < bestFitness
            bestFitness = fitness(1);
            bestPosition = positions(1, :);
        end

        %% Store Convergence
        convergenceCurve(t) = bestFitness;

    end

end
