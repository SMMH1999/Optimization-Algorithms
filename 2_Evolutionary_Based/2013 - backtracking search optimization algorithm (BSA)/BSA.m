function [bestFitness, bestPosition, convergenceCurve] = BSA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Backtracking Search Optimization Algorithm (BSA)
    % -------------------------------------------------------------------------
    % Reference:
    % P. Civicioglu, "Backtracking Search Optimization Algorithm for numerical
    % optimization problems", Applied Mathematics and Computation, 219,
    % 8121–8144, 2013.
    %
    % Algorithm Name: Backtracking Search Optimization Algorithm (BSA)
    % Abbreviation  : BSA
    %
    % Description:
    % BSA is a population-based evolutionary optimization algorithm that
    % utilizes a historical population (swarm-memory) and a backtracking
    % mutation–crossover strategy to explore and exploit the search space.
    % The algorithm is designed for continuous numerical optimization problems.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound of decision variables (scalar or 1×dim vector)
    %   ub        - Upper bound of decision variables (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of individuals)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %               Usage: fitness = objFun(position)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best decision vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (maxItr×1)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   DIM_RATE - Dimensional crossover rate (fixed to 1 as in original code)
    %   F        - Scale factor generated using standard Brownian walk (3*randn)
    %
    % Notes:
    %   - Minimization problem.
    %   - No function evaluation (FEs) counter is used (externally controlled).
    %   - Boundary control follows original stochastic correction strategy.
    % =========================================================================

    %% --------------------------- Initialization -----------------------------

    DIM_RATE = 1;  % As commonly used in BSA

    % Bound handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    % Generate initial population (Eq.1)
    pop = rand(nPop, dim) .* (ub - lb) + lb;

    % Evaluate initial fitness
    fitnessPop = zeros(nPop,1);
    for i = 1:nPop
        fitnessPop(i) = objFun(pop(i,:));
    end

    % Historical population (Eq.2)
    historicalPop = rand(nPop, dim) .* (ub - lb) + lb;

    % Global best initialization
    [bestFitness, idx] = min(fitnessPop);
    bestPosition = pop(idx,:);

    % Convergence curve
    convergenceCurve = zeros(maxItr,1);

    %% ------------------------------ Main Loop -------------------------------
    for t = 1:maxItr

        %% ------------------------- Selection-I ------------------------------

        if rand < rand
            historicalPop = pop;
        end

        historicalPop = historicalPop(randperm(nPop), :);

        % Scale factor (Eq.5) - Brownian walk
        F = 3 * randn;

        % Map matrix (Algorithm-2)
        map = zeros(nPop, dim);

        if rand < rand
            for i = 1:nPop
                u = randperm(dim);
                nChange = ceil(DIM_RATE * rand * dim);
                map(i, u(1:nChange)) = 1;
            end
        else
            for i = 1:nPop
                map(i, randi(dim)) = 1;
            end
        end

        %% -------------------- Recombination (Mutation+Crossover) -----------

        offsprings = pop + (map .* F) .* (historicalPop - pop);

        % Boundary Control (Algorithm-3)
        for i = 1:nPop
            for j = 1:dim
                k = rand < rand;
                if offsprings(i,j) < lb(j)
                    if k
                        offsprings(i,j) = lb(j);
                    else
                        offsprings(i,j) = rand*(ub(j)-lb(j)) + lb(j);
                    end
                end
                if offsprings(i,j) > ub(j)
                    if k
                        offsprings(i,j) = ub(j);
                    else
                        offsprings(i,j) = rand*(ub(j)-lb(j)) + lb(j);
                    end
                end
            end
        end

        %% -------------------------- Selection-II ----------------------------

        fitnessOff = zeros(nPop,1);
        for i = 1:nPop
            fitnessOff(i) = objFun(offsprings(i,:));
        end

        improved = fitnessOff < fitnessPop;
        fitnessPop(improved) = fitnessOff(improved);
        pop(improved,:) = offsprings(improved,:);

        %% -------------------------- Best Update ------------------------------

        [currentBest, idx] = min(fitnessPop);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = pop(idx,:);
        end

        %% ----------------------- Convergence Curve --------------------------

        convergenceCurve(t) = bestFitness;

    end

end
