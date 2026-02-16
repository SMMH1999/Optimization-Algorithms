function [bestFitness, bestPosition, convergenceCurve] = ZOA(lb, ub, dim, nPop, maxItr, objFun)
    % Zebra Optimization Algorithm (ZOA)
    % Abbreviation: ZOA
    %
    % Developed by: Eva Trojovská, Mohammad Dehghani, and Pavel Trojovský
    %
    % Description:
    %   Zebra Optimization Algorithm is a nature-inspired metaheuristic
    %   optimization algorithm simulating zebra foraging and defense strategies
    %   against predators.
    %
    % Inputs:
    %   lb          - scalar or 1xdim vector, lower bound of variables
    %   ub          - scalar or 1xdim vector, upper bound of variables
    %   dim         - integer, number of decision variables (dimensions)
    %   nPop        - integer, number of search agents (population size)
    %   maxItr      - integer, maximum number of iterations
    %   objFun      - function handle, objective function to minimize
    %
    % Outputs:
    %   bestFitness       - best objective function value found
    %   bestPosition      - position vector corresponding to bestFitness
    %   convergenceCurve  - vector of bestFitness at each iteration
    %
    % Tunable Parameters:
    %   None (the algorithm uses default randomization internally)

    %% Boundary setup
    lb = ones(1, dim) .* lb;
    ub = ones(1, dim) .* ub;

    %% INITIALIZATION
    X = zeros(nPop, dim);
    for i = 1:dim
        X(:, i) = lb(i) + rand(nPop, 1) .* (ub(i) - lb(i));
    end

    fit = zeros(nPop, 1);
    for i = 1:nPop
        fit(i) = objFun(X(i, :));
    end

    %% MAIN LOOP
    best_so_far = zeros(maxItr, 1); % Convergence tracker

    for t = 1:maxItr
        %% Update global best
        [currentBest, loc] = min(fit);
        if t == 1
            bestPosition = X(loc, :);
            bestFitness = currentBest;
        elseif currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(loc, :);
        end

        %% PHASE 1: Foraging Behavior
        for i = 1:nPop
            I = round(1 + rand);
            X_new = X(i, :) + rand(1, dim) .* (bestPosition - I .* X(i, :));
            X_new = max(min(X_new, ub), lb); % enforce bounds

            f_new = objFun(X_new);
            if f_new <= fit(i)
                X(i, :) = X_new;
                fit(i) = f_new;
            end
        end

        %% PHASE 2: Defense Strategies Against Predators
        Ps = rand;
        k = randperm(nPop, 1);
        AZ = X(k, :); % attacked zebra

        for i = 1:nPop
            if Ps < 0.5
                % S1: Escape strategy
                R = 0.1;
                X_new = X(i, :) + R * (2 * rand(1, dim) - 1) * (1 - t / maxItr) .* X(i, :);
            else
                % S2: Offensive strategy
                I = round(1 + rand);
                X_new = X(i, :) + rand(1, dim) .* (AZ - I .* X(i, :));
            end

            X_new = max(min(X_new, ub), lb); % enforce bounds
            f_new = objFun(X_new);

            if f_new <= fit(i)
                X(i, :) = X_new;
                fit(i) = f_new;
            end
        end

        %% Update convergence curve
        best_so_far(t) = bestFitness;
    end

    convergenceCurve = best_so_far;
end
