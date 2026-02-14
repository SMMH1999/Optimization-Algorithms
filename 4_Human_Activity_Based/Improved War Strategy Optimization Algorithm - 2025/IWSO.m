function [bestFitness, bestPosition, convergenceCurve] = IWSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Improved War Strategy Optimization (IWSO)
    % =========================================================================
    % Author: Aydilek, İ.B.; Uslu, A.; Kına, C.
    % Ref: Improved War Strategy Optimization with Extreme Learning Machine
    %      for Health Data Classification. Applied Sciences, 2025.
    %
    % Description:
    % IWSO is a population-based metaheuristic inspired by war strategies.
    % It enhances the original War Strategy Optimization (WSO) algorithm by
    % incorporating a Random Opposition-Based Learning (ROBL) mechanism to
    % improve exploration and convergence performance.
    %
    % Input Parameters:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of soldiers)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Output:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best solution vector (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % Algorithm Parameters:
    %   R   - Random control parameter (default = 0.1)
    %   W1  - Adaptive weight vector (initialized to 2)
    %   Wg  - Improvement counter per individual
    %
    % Notes:
    %   - Minimization problem.
    %   - Boundary constraints are strictly enforced.
    %   - Supports scalar or vector bounds.
    % =========================================================================

    %% --------------------- Boundary Handling -------------------------------
    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    %% --------------------- Initialization ----------------------------------
    Positions = rand(nPop, dim) .* (ub - lb) + lb;
    Positions_new = zeros(nPop, dim);

    fitness_old = inf(1, nPop);
    fitness_new = inf(1, nPop);

    bestPosition = zeros(1, dim);
    bestFitness = inf;

    W1 = 2 * ones(1, nPop);
    Wg = zeros(1, nPop);
    R = 0.1;

    convergenceCurve = zeros(1, maxItr);

    %% --------------------- Initial Fitness Evaluation ----------------------
    for i = 1:nPop
        fitness = objFun(Positions(i, :));
        fitness_old(i) = fitness;

        if fitness < bestFitness
            bestFitness = fitness;
            bestPosition = Positions(i, :);
        end
    end

    convergenceCurve(1) = bestFitness;

    %% --------------------- Main Loop ---------------------------------------
    for t = 2:maxItr

        [~, sortedIdx] = sort(fitness_old);
        Co = Positions(sortedIdx(2), :);
        com = randperm(nPop);

        for i = 1:nPop

            RR = rand;

            if RR < R
                D_V = 2 * RR * (bestPosition - Positions(com(i), :)) + ...
                    W1(i) * rand * (Co - Positions(i, :));
            else
                D_V = 2 * RR * (Co - bestPosition) + ...
                    rand * (W1(i) * bestPosition - Positions(i, :));
            end

            Positions_new(i, :) = Positions(i, :) + D_V;

            % Boundary control
            Positions_new(i, :) = max(Positions_new(i, :), lb);
            Positions_new(i, :) = min(Positions_new(i, :), ub);

            fitness = objFun(Positions_new(i, :));
            fitness_new(i) = fitness;

            %% -------- Random Opposition-Based Learning (ROBL) --------------
            xbar_OBL = ub + lb - Positions_new(i, :);
            x_proposed_ROBL = rand .* (xbar_OBL - Positions_new(i, :)) + Positions_new(i, :);

            x_proposed_ROBL = max(x_proposed_ROBL, lb);
            x_proposed_ROBL = min(x_proposed_ROBL, ub);

            fitness_temp = objFun(x_proposed_ROBL);

            if fitness_temp < fitness
                Positions_new(i, :) = x_proposed_ROBL;
                fitness_new(i) = fitness_temp;
                fitness = fitness_temp;
            end

            %% -------- Global Best Update ------------------------------------
            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = Positions_new(i, :);
            end

            %% -------- Greedy Selection --------------------------------------
            if fitness < fitness_old(i)
                Positions(i, :) = Positions_new(i, :);
                fitness_old(i) = fitness;
                Wg(i) = Wg(i) + 1;
                W1(i) = W1(i) * (1 - Wg(i) / maxItr)^2;
            end

        end

        %% -------- Worst Replacement Strategy --------------------------------
        if t < 1000
            [~, worstIdx] = max(fitness_old);
            Positions(worstIdx, :) = rand(1, dim) .* (ub - lb) + lb;
            fitness_old(worstIdx) = objFun(Positions(worstIdx, :));

            if fitness_old(worstIdx) < bestFitness
                bestFitness = fitness_old(worstIdx);
                bestPosition = Positions(worstIdx, :);
            end
        end

        %% -------- Convergence Curve -----------------------------------------
        convergenceCurve(t) = bestFitness;

    end

end
