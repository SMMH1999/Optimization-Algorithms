function [bestFitness, bestPosition, convergenceCurve] = GOA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Grasshopper Optimization Algorithm (GOA)
    % =========================================================================
    % Author: (Original authors of GOA paper)
    % Refactored and formatted for benchmark framework (OAF Project)
    %
    % Description:
    % The Grasshopper Optimization Algorithm (GOA) is a population-based
    % metaheuristic optimization algorithm inspired by the social interaction
    % behavior of grasshoppers. The algorithm consists of two main phases:
    % exploration (global search) and exploitation (local search).
    %
    % =========================================================================
    % Inputs:
    %   lb        - Lower bound of variables (scalar or 1×dim vector)
    %   ub        - Upper bound of variables (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of search agents)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best decision variable vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % =========================================================================
    % Notes:
    % - Minimization problem is assumed.
    % - Boundary constraints are strictly enforced.
    % - No function evaluation (FEs) structure is used (controlled externally).
    % =========================================================================

    %% ========================== Initialization ==============================
    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    % Initialize population
    X = zeros(nPop, dim);
    for i = 1:dim
        X(:, i) = lb(i) + rand(nPop, 1) .* (ub(i) - lb(i));
    end

    % Evaluate initial population
    fit = zeros(nPop, 1);
    for i = 1:nPop
        fit(i) = objFun(X(i, :));
    end

    % Initialize global best
    [bestFitness, bestIdx] = min(fit);
    bestPosition = X(bestIdx, :);

    convergenceCurve = zeros(1, maxItr);

    %% ============================ Main Loop ================================
    for t = 1:maxItr

        % -------- Update Global Best --------
        [currentBest, currentIdx] = min(fit);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(currentIdx, :);
        end

        %% ================= Phase 1: Exploration ============================
        for i = 1:nPop

            if rand < 0.5
                I = round(1 + rand(1,1));
                RAND = rand(1,1);
            else
                I = round(1 + rand(1,dim));
                RAND = rand(1,dim);
            end

            X_P1 = X(i,:) + RAND .* (bestPosition - I .* X(i,:));   % Eq. (4)

            % Boundary control
            X_P1 = max(X_P1, lb);
            X_P1 = min(X_P1, ub);

            % Fitness evaluation
            F_P1 = objFun(X_P1);

            % Greedy selection (Eq. 5)
            if F_P1 < fit(i)
                X(i,:) = X_P1;
                fit(i) = F_P1;
            end
        end

        %% ================= Phase 2: Exploitation ===========================
        for i = 1:nPop

            X_P2 = X(i,:) + (1 - 2*rand(1,1)) .* ...
                ( lb./t + rand(1,1).*(ub./t - lb./t) );          % Eq. (6)

            % Boundary control (scaled then original bounds)
            X_P2 = max(X_P2, lb./t);
            X_P2 = min(X_P2, ub./t);
            X_P2 = max(X_P2, lb);
            X_P2 = min(X_P2, ub);

            % Fitness evaluation
            F_P2 = objFun(X_P2);

            % Greedy selection (Eq. 7)
            if F_P2 < fit(i)
                X(i,:) = X_P2;
                fit(i) = F_P2;
            end
        end

        %% ================= Convergence Curve ===============================
        convergenceCurve(t) = bestFitness;

    end

end
