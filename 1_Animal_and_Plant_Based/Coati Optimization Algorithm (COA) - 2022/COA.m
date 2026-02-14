function [bestFitness, bestPosition, convergenceCurve] = COA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % COA - Coati Optimization Algorithm
    % =========================================================================
    % Designed by: Zeinab Montazeri
    % Developed by: Mohammad Dehghani and Pavel Trojovský
    %
    % Refactored to benchmark single-function format.
    %
    % -------------------------------------------------------------------------
    % Description:
    % COA is a population-based metaheuristic inspired by the hunting behavior
    % of coatis. The algorithm consists of two main phases:
    %   1) Exploration Phase: Hunting and attacking strategy on iguana.
    %   2) Exploitation Phase: Escaping from predators.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of search agents)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best solution vector found
    %   convergenceCurve - Best fitness value at each iteration
    %
    % -------------------------------------------------------------------------
    % Notes:
    % - Minimization problem is assumed.
    % - Boundary constraints are strictly enforced.
    % - Algorithm logic preserved exactly from the original implementation.
    % =========================================================================

    %% --------------------------- Bound Handling ----------------------------
    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    %% --------------------------- Initialization ----------------------------
    X = zeros(nPop, dim);
    for i = 1:dim
        X(:, i) = lb(i) + rand(nPop, 1) .* (ub(i) - lb(i));
    end

    fit = zeros(nPop, 1);
    for i = 1:nPop
        fit(i) = objFun(X(i, :));
    end

    bestFitness = inf;
    bestPosition = zeros(1, dim);
    convergenceCurve = zeros(1, maxItr);

    %% ============================ Main Loop ================================
    for t = 1:maxItr

        %% --------------------- Update Global Best --------------------------
        [currentBest, location] = min(fit);
        if t == 1
            bestFitness = currentBest;
            bestPosition = X(location, :);
        elseif currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(location, :);
        end

        %% ===================================================================
        %% Phase 1: Hunting and Attacking Strategy (Exploration Phase)
        %% ===================================================================

        % First half of population
        for i = 1:nPop/2

            iguana = bestPosition;
            I = round(1 + rand);

            X_P1 = X(i, :) + rand .* (iguana - I .* X(i, :));   % Eq. (4)
            X_P1 = max(X_P1, lb);
            X_P1 = min(X_P1, ub);

            F_P1 = objFun(X_P1);                                 % Eq. (7)

            if F_P1 < fit(i)
                X(i, :) = X_P1;
                fit(i) = F_P1;
            end
        end

        % Second half of population
        for i = (1 + nPop/2):nPop

            iguana = lb + rand .* (ub - lb);                     % Eq. (5)
            F_HL = objFun(iguana);
            I = round(1 + rand);

            if fit(i) > F_HL
                X_P1 = X(i, :) + rand .* (iguana - I .* X(i, :)); % Eq. (6)
            else
                X_P1 = X(i, :) + rand .* (X(i, :) - iguana);      % Eq. (6)
            end

            X_P1 = max(X_P1, lb);
            X_P1 = min(X_P1, ub);

            F_P1 = objFun(X_P1);                                 % Eq. (7)

            if F_P1 < fit(i)
                X(i, :) = X_P1;
                fit(i) = F_P1;
            end
        end

        %% ===================================================================
        %% Phase 2: Escaping from Predators (Exploitation Phase)
        %% ===================================================================

        for i = 1:nPop

            LO_LOCAL = lb / t;                                    % Eq. (9)
            HI_LOCAL = ub / t;                                    % Eq. (10)

            X_P2 = X(i, :) + (1 - 2 * rand) .* ...
                (LO_LOCAL + rand .* (HI_LOCAL - LO_LOCAL));     % Eq. (8)

            X_P2 = max(X_P2, LO_LOCAL);
            X_P2 = min(X_P2, HI_LOCAL);

            F_P2 = objFun(X_P2);                                  % Eq. (11)

            if F_P2 < fit(i)
                X(i, :) = X_P2;
                fit(i) = F_P2;
            end
        end

        %% ------------------------ Convergence ------------------------------
        convergenceCurve(t) = bestFitness;

    end

end
