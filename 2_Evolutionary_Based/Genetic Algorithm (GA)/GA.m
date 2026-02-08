function [bestFitness, bestPosition, convergenceCurve] = GA(LB, UB, Dim, popSize, maxItr, Cost_Function)
    % GA - Simple Genetic Algorithm (real-coded, rank selection, 1-point crossover, uniform mutation)

    %% Parameters
    pc       = 0.8;    % crossover probability
    pm       = 0.01;   % mutation probability per gene
    eliteCnt = 1;      % number of elites

    % Expand bounds if scalar
    if isscalar(LB), LB = LB * ones(1, Dim); end
    if isscalar(UB), UB = UB * ones(1, Dim); end

    convergenceCurve = inf(1, maxItr);
    fitness = inf(popSize, 1);

    %% Initial population
    pop = rand(popSize, Dim) .* (UB - LB) + LB;

    %% Initial evaluation
    for i = 1:popSize
        fitness(i) = Cost_Function(pop(i,:));
    end

    [fitness, order] = sort(fitness);
    pop = pop(order, :);

    bestFitness  = fitness(1);
    bestPosition = pop(1,:);

    %% Main loop
    for t = 1:maxItr

        convergenceCurve(t) = bestFitness;

        %% Elitism
        elites = pop(1:eliteCnt, :);

        %% Selection (roulette on rank)
        ranks   = popSize - (1:popSize)' + 1;   % best rank = popSize
        selProb = ranks / sum(ranks);
        cumProb = cumsum(selProb);

        matingPool = zeros(popSize - eliteCnt, Dim);
        for k = 1:size(matingPool,1)
            r = rand;
            idx = find(cumProb >= r, 1, 'first');
            matingPool(k,:) = pop(idx,:);
        end

        %% Crossover (1-point)
        offspring = matingPool;
        for i = 1:2:size(matingPool,1)-1
            if rand < pc
                cp = randi([1, Dim-1]);
                offspring(i, cp+1:end)   = matingPool(i+1, cp+1:end);
                offspring(i+1, cp+1:end) = matingPool(i,   cp+1:end);
            end
        end

        %% Mutation (uniform)
        for i = 1:size(offspring,1)
            for j = 1:Dim
                if rand < pm
                    offspring(i,j) = LB(j) + rand * (UB(j) - LB(j));
                end
            end
        end

        %% Form next generation
        pop = [elites; offspring];

        %% Boundary control and evaluation
        pop = max(min(pop, UB), LB);

        for i = 1:popSize
            fitness(i) = Cost_Function(pop(i,:));
        end

        [fitness, order] = sort(fitness);
        pop = pop(order, :);

        %% Update global best
        if fitness(1) < bestFitness
            bestFitness  = fitness(1);
            bestPosition = pop(1,:);
        end
    end

    % Enforce monotonic convergence
    convergenceCurve = cummin(convergenceCurve);
end
