function [bestFitness, bestPosition, convergenceCurve] = OOA(lb, ub, dim, nPop, maxItr, objFun)
    %% Optimal Ocean Algorithm (OOA)
    %   Authors: Pavel Trojovský and Mohammad Dehghani
    %   MATLAB Implementation: Refactored Benchmark Version
    %
    % Description:
    %   OOA is a population-based optimization algorithm inspired by
    %   ocean fish movement and hunting strategies.
    %
    % Inputs:
    %   nPop      - number of search agents (integer)
    %   maxItr    - maximum number of iterations (integer)
    %   lb        - lower bound (scalar or 1xdim vector)
    %   ub        - upper bound (scalar or 1xdim vector)
    %   dim       - number of decision variables (integer)
    %   objFun    - handle to objective function, e.g., @(x) sum(x.^2)
    %
    % Outputs:
    %   bestFitness      - best objective function value found (scalar)
    %   bestPosition     - position vector corresponding to bestFitness (1xdim)
    %   convergenceCurve - array of bestFitness at each iteration (1xmaxItr)
    %
    % Tunable Parameters:
    %   nPop   - number of ocean fish (population size)
    %   maxItr - maximum iterations
    %
    % Usage Example:
    %   [bestF, bestX, curve] = OOA(30, 500, -10, 10, 5, @(x) sum(x.^2));

    %% Ensure bounds are vectors
    lb = ones(1, dim) .* lb;
    ub = ones(1, dim) .* ub;

    %% ---------------------------- Initialization ----------------------------
    X = zeros(nPop, dim);       % Population positions
    fit = zeros(nPop, 1);       % Fitness values

    for i = 1:dim
        X(:, i) = lb(i) + rand(nPop, 1) .* (ub(i) - lb(i));
    end

    for i = 1:nPop
        fit(i) = objFun(X(i, :));
    end

    best_so_far = zeros(1, maxItr);  % Track best fitness per iteration

    %% --------------------------- Main Optimization Loop ---------------------------
    for t = 1:maxItr
        % Find current best solution
        [Fbest, bestIdx] = min(fit);
        if t == 1
            xbest = X(bestIdx, :);
            fbest = Fbest;
        elseif Fbest < fbest
            fbest = Fbest;
            xbest = X(bestIdx, :);
        end

        % Update each search agent
        for i = 1:nPop
            %% Phase 1: Exploration (hunting the fish)
            fishPos = find(fit < fit(i));
            if isempty(fishPos)
                selectedFish = xbest;
            else
                if rand < 0.5
                    selectedFish = xbest;
                else
                    k = randperm(numel(fishPos), 1);
                    selectedFish = X(fishPos(k), :);
                end
            end

            I = round(1 + rand);
            X_new = X(i, :) + rand(1, dim) .* (selectedFish - I .* X(i, :));
            X_new = max(X_new, lb);
            X_new = min(X_new, ub);

            % Update position if improved
            fit_new = objFun(X_new);
            if fit_new < fit(i)
                X(i, :) = X_new;
                fit(i) = fit_new;
            end

            %% Phase 2: Exploitation (carrying fish to suitable positions)
            X_new = X(i, :) + (lb + rand(1, dim) .* (ub - lb)) / t;
            X_new = max(X_new, lb);
            X_new = min(X_new, ub);

            fit_new = objFun(X_new);
            if fit_new < fit(i)
                X(i, :) = X_new;
                fit(i) = fit_new;
            end
        end

        % Record convergence
        best_so_far(t) = fbest;
    end

    %% --------------------------- Return Results ---------------------------
    bestFitness = fbest;
    bestPosition = xbest;
    convergenceCurve = best_so_far;

end
