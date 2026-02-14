function [Best_Fitness, Best_Position, Convergence_Curve] = DE(LB, UB, Dim, Pop_size, Max_iter, Cost_Function)

    %% Expand bounds if scalar
    if isscalar(UB)
        UB = repmat(UB, 1, Dim);
        LB = repmat(LB, 1, Dim);
    end

    %% Initialization
    Population = Population_Generator(Pop_size, Dim, UB, LB);
    Fitness = zeros(1, Pop_size);
    Convergence_Curve = inf(1, Max_iter);

    % Evaluate initial population
    for i = 1:Pop_size
        Fitness(i) = Cost_Function(Population(i, :));
    end

    %% DE Parameters
    F  = 0.5;   % Mutation factor
    CR = 0.9;   % Crossover rate

    %% Main loop
    for itr = 1:Max_iter

        for i = 1:Pop_size

            %% Mutation (DE/rand/1)
            idxs = randperm(Pop_size, 3);
            while any(idxs == i)
                idxs = randperm(Pop_size, 3);
            end

            x1 = Population(idxs(1), :);
            x2 = Population(idxs(2), :);
            x3 = Population(idxs(3), :);

            mutant = x1 + F * (x2 - x3);

            % Boundary control
            mutant = max(min(mutant, UB), LB);

            %% Crossover (binomial)
            trial = Population(i, :);
            j_rand = randi(Dim);

            for j = 1:Dim
                if rand <= CR || j == j_rand
                    trial(j) = mutant(j);
                end
            end

            %% Selection
            trial_fitness = Cost_Function(trial);

            if trial_fitness < Fitness(i)
                Population(i, :) = trial;
                Fitness(i) = trial_fitness;
            end
        end

        %% Store best
        [Best_Fitness, best_idx] = min(Fitness);
        Convergence_Curve(itr) = Best_Fitness;
    end

    %% Final best position
    Best_Position = Population(best_idx, :);
end
