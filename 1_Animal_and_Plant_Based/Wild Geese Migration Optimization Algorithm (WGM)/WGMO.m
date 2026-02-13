function [bestFitness, bestPosition, convergenceCurve] = WGMO(lb, ub, dim, nPop, maxItr, objFun)
    %==========================================================================
    % Wild Geese Migration Optimization (WGMO)
    %--------------------------------------------------------------------------
    % Author: Your Name (if available)
    % Description: Implementation of the Wild Geese Migration Optimization
    %              algorithm for continuous optimization problems.
    %--------------------------------------------------------------------------
    % INPUTS:
    %   lb       - Lower bound (scalar or 1xdim vector)
    %   ub       - Upper bound (scalar or 1xdim vector)
    %   dim      - Number of dimensions (variables)
    %   nPop     - Population size (number of geese)
    %   maxItr   - Maximum number of iterations
    %   objFun   - Handle to objective function: fitness = objFun(position)
    %--------------------------------------------------------------------------
    % OUTPUTS:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Position of the best solution
    %   convergenceCurve - Array containing the best fitness at each iteration
    %--------------------------------------------------------------------------
    % TUNABLE PARAMETERS:
    %   alpha  - Inertia / scaling factor (default: 0.9)
    %   beta   - Learning rate toward best goose (default: 0.1)
    %   gamma  - Randomization factor (default: 0.1)
    %==========================================================================

    %% Parameters
    alpha = 0.9;  % scaling factor
    beta  = 0.1;  % learning rate
    gamma = 0.1;  % randomization factor

    % Ensure lb and ub are vectors
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% Initialization
    geese = rand(nPop, dim) .* (ub - lb) + lb;  % population positions
    fitness = zeros(nPop, 1);                   % fitness values
    convergenceCurve = zeros(maxItr, 1);        % store best fitness per iteration

    %% Main Optimization Loop
    for t = 1:maxItr
        % Evaluate fitness
        for i = 1:nPop
            fitness(i) = objFun(geese(i, :));
        end

        % Sort geese based on fitness (ascending)
        [fitness, idx] = sort(fitness);
        geese = geese(idx, :);

        % Update positions of all geese
        best_goose = geese(1, :);  % best goose in current iteration
        for i = 1:nPop
            geese(i, :) = alpha * geese(i, :) ...
                + beta * rand(1, dim) .* (best_goose - geese(i, :)) ...
                + gamma * rand(1, dim) .* (ub - lb);
            % Boundary check
            geese(i, :) = min(max(geese(i, :), lb), ub);
        end

        % Update convergence curve
        convergenceCurve(t) = fitness(1);
    end

    %% Output best solution
    bestFitness = fitness(1);
    bestPosition = geese(1, :);

end
