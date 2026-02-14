function [bestSolution, bestFitness] = GOA(dim, N, maxIter, lb, ub, objFun)
    % Goat Optimization Algorithm (GOA) for Single-Objective Optimization
    %
    % Inputs:
    % dim      - Problem dimension (scalar)
    % N        - Population size (number of goats)
    % maxIter  - Maximum number of iterations
    % lb       - Lower bound (vector or scalar)
    % ub       - Upper bound (vector or scalar)
    % objFun   - Handle to the objective function: f = objFun(x)
    %
    % Outputs:
    % bestSolution - Best solution found (1 x dim)
    % bestFitness  - Best fitness value

    % Initialize population
    pop = repmat(lb, N, 1) + rand(N, dim) .* (repmat(ub - lb, N, 1));
    fitness = arrayfun(@(i) objFun(pop(i, :)), 1:N)';
    [bestFitness, bestIdx] = min(fitness);
    bestSolution = pop(bestIdx, :);

    % Coefficients
    alpha = 0.05; % Exploration
    beta  = 0.5;  % Exploitation
    J     = 0.1;  % Jump probability

    % Main loop
    for t = 1:maxIter
        for i = 1:N
            r1 = rand(); r2 = rand(); r3 = rand();

            % Exploration
            if r1 < 0.5
                pop(i, :) = pop(i, :) + alpha * randn(1, dim) .* (ub - lb);
            end

            % Exploitation
            if r2 >= 0.5
                pop(i, :) = pop(i, :) + beta * (bestSolution - pop(i, :));
            end

            % Jump strategy
            if r3 < J
                randIdx = randi([1, N]);
                pop(i, :) = pop(i, :) + J * (pop(randIdx, :) - pop(i, :));
            end

            % Apply bounds
            pop(i, :) = min(max(pop(i, :), lb), ub);
        end

        % Evaluate fitness
        fitness = arrayfun(@(i) objFun(pop(i, :)), 1:N)';

        % Parasite avoidance (replace weakest 20%)
        [~, sortedIdx] = sort(fitness);
        weakest = sortedIdx(end - floor(0.2*N) + 1:end);
        pop(weakest, :) = repmat(lb, length(weakest), 1) + rand(length(weakest), dim) .* (repmat(ub - lb, length(weakest), 1));

        % Update best solution
        [currentBestFitness, bestIdx] = min(fitness);
        if currentBestFitness < bestFitness
            bestFitness = currentBestFitness;
            bestSolution = pop(bestIdx, :);
        end
    end
end
