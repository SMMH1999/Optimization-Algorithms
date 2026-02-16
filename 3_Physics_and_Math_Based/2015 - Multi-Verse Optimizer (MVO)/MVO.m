function [bestFitness, bestPosition, convergenceCurve] = MVO(lb, ub, dim, nPop, maxItr, objFun)
    %% Multi-Verse Optimizer (MVO) - Unified Benchmark Version
    %  Abbreviation: MVO
    %
    %  Author: Seyedali Mirjalili
    %  Reference Paper:
    %    S. Mirjalili, S. M. Mirjalili, A. Hatamlou,
    %    Multi-Verse Optimizer: a nature-inspired algorithm for global optimization,
    %    Neural Computing and Applications, 2015,
    %    DOI: http://dx.doi.org/10.1007/s00521-015-1870-7
    %
    %  Description:
    %    Multi-Verse Optimizer (MVO) is a population-based metaheuristic algorithm
    %    inspired by multiverse concepts and wormholes.
    %
    %  Inputs:
    %    lb       : Lower bound(s) [scalar or 1 x dim vector]
    %    ub       : Upper bound(s) [scalar or 1 x dim vector]
    %    dim      : Number of decision variables
    %    nPop     : Population size
    %    maxItr   : Maximum number of iterations
    %    objFun   : Function handle to the objective function
    %
    %  Outputs:
    %    bestFitness       : Best fitness value found
    %    bestPosition      : Position of the best solution
    %    convergenceCurve  : Best fitness at each iteration

    %% Initialization
    bestPosition = zeros(1, dim);
    bestFitness = inf;
    population = initializePopulation(nPop, dim, lb, ub);

    WEP_Max = 1;
    WEP_Min = 0.2;
    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for iter = 1:maxItr
        WEP = WEP_Min + iter * ((WEP_Max - WEP_Min) / maxItr);
        TDR = 1 - (iter^(1/6) / maxItr^(1/6));

        fitness = zeros(nPop,1);

        % Evaluate fitness and update best
        for i = 1:nPop
            population(i,:) = max(min(population(i,:), ub), lb); % boundary check
            fitness(i) = objFun(population(i,:));
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = population(i,:);
            end
        end

        % Sort universes
        [sortedFitness, idx] = sort(fitness);
        sortedPopulation = population(idx,:);
        normFitness = normr(sortedFitness(:));

        % Elite preservation
        population(1,:) = sortedPopulation(1,:);

        % Update positions
        for i = 2:nPop
            for j = 1:dim
                if rand < normFitness(i)
                    wh_idx = RouletteWheelSelection(-sortedFitness);
                    if wh_idx == -1
                        wh_idx = 1;
                    end
                    population(i,j) = sortedPopulation(wh_idx,j);
                end

                if rand < WEP
                    r = rand;
                    if isscalar(lb)
                        delta = (ub - lb) * rand + lb;
                    else
                        delta = (ub(j) - lb(j)) * rand + lb(j);
                    end
                    if r < 0.5
                        population(i,j) = bestPosition(j) + TDR * delta;
                    else
                        population(i,j) = bestPosition(j) - TDR * delta;
                    end
                end
            end
        end

        convergenceCurve(iter) = bestFitness;
    end

end

%% --- Helper Functions ---
function pop = initializePopulation(nPop, dim, ub, lb)
    if isscalar(ub)
        pop = rand(nPop, dim) * (ub - lb) + lb;
    else
        pop = zeros(nPop, dim);
        for i = 1:dim
            pop(:,i) = rand(nPop,1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end

function idx = RouletteWheelSelection(weights)
    acc = cumsum(weights);
    p = rand() * acc(end);
    idx = -1;
    for i = 1:length(acc)
        if acc(i) > p
            idx = i;
            break;
        end
    end
end
