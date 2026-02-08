function [bestScore, bestPos, curve] = GA(LB, UB, Dim, popSize, maxItr, objective)
    % GA_RS: Genetic Algorithm with interface matching RS
    % Inputs:
    %   LB       - Lower bounds (scalar or vector)
    %   UB       - Upper bounds (scalar or vector)
    %   Dim      - Number of dimensions
    %   popSize  - Population size
    %   maxItr   - Maximum number of iterations
    %   objective - Function handle: f = objective(x), where x is a row vector
    %
    % Outputs:
    %   bestScore - Best fitness found
    %   bestPos   - Best solution (position)
    %   curve     - Convergence curve (best fitness at each iteration)
    
    % --- Normalize bounds ---
    if isscalar(LB), LB = LB * ones(1, Dim); end
    if isscalar(UB), UB = UB * ones(1, Dim); end
    
    %% -------------------- Initialize --------------------
    bestScore = inf;
    bestPos   = zeros(1, Dim);
    curve     = zeros(maxItr, 1);
    
    % Initial population
    pop = rand(popSize, Dim) .* (UB - LB) + LB;
    
    % Evaluate initial population
    fitness = zeros(popSize, 1);
    for i = 1:popSize
        fitness(i) = objective(pop(i, :));
    end
    
    % --- Sort and update best ---
    [fitness, idx] = sort(fitness);
    pop = pop(idx, :);
    bestScore = fitness(1);
    bestPos   = pop(1, :);
    
    % --- GA Parameters ---
    pc = 0.8;        % Crossover probability
    pm = 0.01;       % Mutation probability per gene
    eliteCnt = 1;    % Elitism: keep best individual
    
    %% -------------------- Main loop --------------------
    for t = 1:maxItr
        curve(t) = bestScore;
        
        % --- Elitism: keep best individual ---
        elites = pop(1:eliteCnt, :);
        
        % --- Rank-based selection (Roulette wheel on ranks) ---
        ranks = popSize - (1:popSize)' + 1;  % Best has highest rank
        selProb = ranks / sum(ranks);
        cumProb = cumsum(selProb);
        matingPool = zeros(popSize - eliteCnt, Dim);
        for k = 1:size(matingPool, 1)
            r = rand;
            idx = find(cumProb >= r, 1, 'first');
            matingPool(k, :) = pop(idx, :);
        end
        
        % --- One-point crossover ---
        offspring = matingPool;
        for i = 1:2:size(matingPool, 1) - 1
            if rand < pc
                cp = randi([1, Dim - 1]);
                offspring(i, cp+1:end) = matingPool(i+1, cp+1:end);
                offspring(i+1, cp+1:end) = matingPool(i, cp+1:end);
            end
        end
        
        % --- Uniform mutation ---
        for i = 1:size(offspring, 1)
            for j = 1:Dim
                if rand < pm
                    offspring(i, j) = LB(j) + rand * (UB(j) - LB(j));
                end
            end
        end


        wSize = 100;

        bFit = min(fitness);
        wFit = max(fitness);

        % normalFit = 1 - ((fitness - bFit) ./ (wFit - bFit));

        normalPos = (pop - LB) ./ (UB - LB);
        normalPos = normalPos * 2 - 1;

        normalFit = 1 - (fitness - bFit) ./ (wFit - bFit); 


        normalFitPos = normalPos .* normalFit;
        normalFitPos =+ normalFitPos/wSize;
        

        % for i = 1 : popSize
        %     hitsogramPop(it) = {pop .* normalFit};
        % end


        if mod(t, wSize) == 0

            % objective(pop(1, :) - mean(normalFitPos(1:5, :)))
            % objective(pop(1, :))
            % objective(pop(1, :) + mean(normalFitPos(1:5, :)))
            popT = normalFitPos .* UB;
            fitT = objective((popT + pop)/2)';
            [fitT, idx] = sort(fitT);
            popT = popT(idx, :);
            disp('1');

        end


        % --- Form new population: elites + offspring ---
        pop = [elites; offspring];
        
        % --- Clip to bounds ---
        pop = max(min(pop, UB), LB);
        
        % --- Evaluate new population ---
        fitness = zeros(popSize, 1);
        for i = 1:popSize
            fitness(i) = objective(pop(i, :));
        end
        
        % --- Sort and update best ---
        [fitness, idx] = sort(fitness);
        pop = pop(idx, :);
        
        % --- Update global best ---
        if fitness(1) < bestScore
            bestScore = fitness(1);
            bestPos   = pop(1, :);
        end
    end
    
    % --- Final convergence curve (cumulative minimum) ---
    curve = cummin(curve);
end