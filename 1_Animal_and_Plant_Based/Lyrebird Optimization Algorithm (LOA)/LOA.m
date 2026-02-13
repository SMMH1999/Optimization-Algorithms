function [bestFitness, bestPosition, convergenceCurve] = LOA(lb, ub, dim, nPop, maxItr, objFun)
    %_______________________________________________________________________________________
    % Lyrebird Optimization Algorithm (LOA)
    %
    % Developed in MATLAB
    %
    % Description:
    %   LOA is a population-based metaheuristic optimization algorithm inspired by
    %   the singing and mating behavior of lyrebirds.
    %
    % Inputs:
    %   lb        - Lower bound (scalar or vector) [1 x dim]
    %   ub        - Upper bound (scalar or vector) [1 x dim]
    %   dim       - Number of decision variables
    %   nPop      - Population size (number of lyrebirds)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Handle to the objective function, e.g., @(x) sphere(x)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Position vector corresponding to best fitness
    %   convergenceCurve - Best fitness value at each iteration
    %
    % Tunable Parameters:
    %   crossover_rate - Probability of performing crossover (default 0.8)
    %   mutation_rate  - Probability of performing mutation (default 0.1)
    %   sigma          - Mutation standard deviation (default 0.1)
    %_______________________________________________________________________________________

    %% Parameters
    crossover_rate = 0.8;   % Crossover probability
    mutation_rate  = 0.1;   % Mutation probability
    sigma          = 0.1;   % Mutation standard deviation

    %% Initialization
    % Handle scalar or vector bounds
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    % Initialize population
    lyrebirds = rand(nPop, dim) .* (ub - lb) + lb;  % Random positions
    fitness = zeros(nPop, 1);                        % Fitness values

    % Convergence history
    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for iteration = 1:maxItr
        % Evaluate fitness
        for i = 1:nPop
            fitness(i) = objFun(lyrebirds(i, :));
        end

        % Sort lyrebirds based on fitness (minimization)
        [fitness, sortedIdx] = sort(fitness);
        lyrebirds = lyrebirds(sortedIdx, :);

        % Update global best
        bestPosition = lyrebirds(1, :);
        bestFitness  = fitness(1);
        convergenceCurve(iteration) = bestFitness;

        % Selection: keep top half
        lyrebirds = lyrebirds(1:round(nPop/2), :);

        %% Crossover
        num_crossovers = round(crossover_rate * size(lyrebirds,1));
        for i = 1:num_crossovers
            % Select two parents randomly
            parentIdx = randi([1, size(lyrebirds,1)], 2, 1);
            parent1 = lyrebirds(parentIdx(1), :);
            parent2 = lyrebirds(parentIdx(2), :);

            % Blend crossover
            alpha = rand();
            child1 = alpha * parent1 + (1 - alpha) * parent2;
            child2 = (1 - alpha) * parent1 + alpha * parent2;

            % Add children
            lyrebirds = [lyrebirds; child1; child2];
        end

        %% Mutation
        num_mutations = round(mutation_rate * nPop);
        for i = 1:num_mutations
            idx = randi([1, size(lyrebirds,1)]);
            mutation = sigma * randn(1, dim);
            lyrebirds(idx, :) = min(max(lyrebirds(idx, :) + mutation, lb), ub);
        end

        % Ensure population size does not exceed nPop (optional trimming)
        if size(lyrebirds,1) > nPop
            lyrebirds = lyrebirds(1:nPop, :);
        end
    end
end
