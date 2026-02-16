function [bestFitness, bestPosition, convergenceCurve] = GOA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Goat Optimization Algorithm (GOA)
    % =========================================================================
    % Author: (Original GOA Author - If Available)
    %
    % Description:
    % Goat Optimization Algorithm (GOA) is a population-based metaheuristic
    % optimizer designed for single-objective optimization problems. The
    % algorithm simulates goat behavioral mechanisms including exploration,
    % exploitation toward the best solution, random jumping behavior, and
    % parasite avoidance (weakest replacement strategy).
    %
    % =========================================================================
    % Inputs:
    % lb        - Lower bound (scalar or 1×dim vector)
    % ub        - Upper bound (scalar or 1×dim vector)
    % dim       - Number of decision variables
    % nPop      - Population size (number of goats)
    % maxItr    - Maximum number of iterations
    % objFun    - Objective function handle (fitness = objFun(position))
    %
    % =========================================================================
    % Outputs:
    % bestFitness      - Best fitness value found (scalar)
    % bestPosition     - Best solution vector (1×dim)
    % convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % =========================================================================
    % Tunable Parameters:
    % alpha  - Exploration coefficient (default = 0.05)
    % beta   - Exploitation coefficient (default = 0.5)
    % J      - Jump probability and jump intensity (default = 0.1)
    %
    % =========================================================================
    % Algorithm Structure:
    % 1. Initialization
    % 2. Fitness Evaluation
    % 3. Main Optimization Loop
    %    - Exploration
    %    - Exploitation
    %    - Jump Strategy
    %    - Boundary Control
    %    - Parasite Avoidance
    %    - Global Best Update
    % 4. Convergence Curve Recording
    % =========================================================================

    %% ------------------------- Initialization ------------------------------

    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    pop = repmat(lb, nPop, 1) + rand(nPop, dim) .* repmat((ub - lb), nPop, 1);

    fitness = zeros(nPop, 1);
    for i = 1:nPop
        fitness(i) = objFun(pop(i, :));
    end

    [bestFitness, bestIdx] = min(fitness);
    bestPosition = pop(bestIdx, :);

    convergenceCurve = zeros(1, maxItr);

    % GOA Parameters
    alpha = 0.05;   % Exploration factor
    beta  = 0.5;    % Exploitation factor
    J     = 0.1;    % Jump probability

    %% ------------------------- Main Loop -----------------------------------

    for t = 1:maxItr

        for i = 1:nPop

            r1 = rand();
            r2 = rand();
            r3 = rand();

            % ---------------- Exploration ----------------
            if r1 < 0.5
                pop(i, :) = pop(i, :) + alpha * randn(1, dim) .* (ub - lb);
            end

            % ---------------- Exploitation ----------------
            if r2 >= 0.5
                pop(i, :) = pop(i, :) + beta * (bestPosition - pop(i, :));
            end

            % ---------------- Jump Strategy ----------------
            if r3 < J
                randIdx = randi([1, nPop]);
                pop(i, :) = pop(i, :) + J * (pop(randIdx, :) - pop(i, :));
            end

            % ---------------- Boundary Control ----------------
            pop(i, :) = min(max(pop(i, :), lb), ub);

        end

        %% ---------------- Fitness Evaluation ----------------

        for i = 1:nPop
            fitness(i) = objFun(pop(i, :));
        end

        %% ---------------- Parasite Avoidance ----------------

        [~, sortedIdx] = sort(fitness);
        nWeak = floor(0.2 * nPop);

        if nWeak > 0
            weakestIdx = sortedIdx(end - nWeak + 1:end);

            pop(weakestIdx, :) = repmat(lb, nWeak, 1) + ...
                rand(nWeak, dim) .* repmat((ub - lb), nWeak, 1);

            for k = 1:nWeak
                fitness(weakestIdx(k)) = objFun(pop(weakestIdx(k), :));
            end
        end

        %% ---------------- Best Update ----------------

        [currentBestFitness, bestIdx] = min(fitness);

        if currentBestFitness < bestFitness
            bestFitness = currentBestFitness;
            bestPosition = pop(bestIdx, :);
        end

        %% ---------------- Convergence Curve ----------------

        convergenceCurve(t) = bestFitness;

    end

end
