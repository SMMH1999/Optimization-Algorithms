function [bestFitness, bestPosition, convergenceCurve] = HGS(lb, ub, dim, nPop, maxItr, objFun)
    %% Hunger Games Search (HGS) - Refactored for Benchmark
    % Author: Yutao Yang, Huiling Chen, Ali Asghar Heidari, Amir H Gandomi
    %
    % Description:
    %   Hunger Games Search (HGS) is a bio-inspired optimization algorithm
    %   for global minimization problems.
    %
    % Inputs:
    %   lb        : lower bound (scalar or 1xdim vector)
    %   ub        : upper bound (scalar or 1xdim vector)
    %   dim       : number of decision variables
    %   nPop      : population size
    %   maxItr    : maximum number of iterations
    %   objFun    : handle to objective function, e.g., @(x) sum(x.^2)
    %
    % Outputs:
    %   bestFitness      : best (minimal) objective function value found
    %   bestPosition     : corresponding decision variable vector
    %   convergenceCurve : bestFitness value at each iteration

    %% Initialization
    bestPosition = zeros(1, dim);
    Destination_fitness = inf;          % Minimization
    Worstest_fitness = -inf;

    X = initializePopulation(nPop, dim, ub, lb);
    convergenceCurve = zeros(1, maxItr);

    hungry = zeros(1, nPop);
    VC1 = ones(nPop, 1);
    weight3 = ones(nPop, dim);
    weight4 = ones(nPop, dim);
    tempPosition = zeros(nPop, dim);

    %% Main Loop
    for it = 1:maxItr
        VC2 = 0.03;   % variation control parameter
        sumHungry = 0;

        % Evaluate fitness and enforce boundaries
        AllFitness = zeros(nPop,1);
        for i = 1:nPop
            X(i,:) = min(max(X(i,:), lb), ub);
            AllFitness(i) = objFun(X(i,:));
        end

        [AllFitnessSorted, IndexSorted] = sort(AllFitness);
        bestIterFitness = AllFitnessSorted(1);
        worstIterFitness = AllFitnessSorted(end);

        % Update global best
        if bestIterFitness < Destination_fitness
            bestPosition = X(IndexSorted(1), :);
            Destination_fitness = bestIterFitness;
            count = 0;
        end

        if worstIterFitness > Worstest_fitness
            Worstest_fitness = worstIterFitness;
        end

        % Update hungry values
        count = 0;
        for i = 1:nPop
            VC1(i) = sech(abs(AllFitness(i) - Destination_fitness));
            if AllFitness(i) == Destination_fitness
                hungry(i) = 0;
                count = count + 1;
                tempPosition(count,:) = X(i,:);
            else
                temprand = rand();
                c = (AllFitness(i) - Destination_fitness) / ...
                    (Worstest_fitness - Destination_fitness) * temprand * 2 .* (ub - lb);
                b = max(c, 100*(1+temprand));
                hungry(i) = hungry(i) + b;
                sumHungry = sumHungry + hungry(i);
            end
        end

        % Compute hungry weights
        for i = 1:nPop
            for j = 1:dim
                weight3(i,j) = (1 - exp(-abs(hungry(i) - sumHungry))) * rand() * 2;
                if rand() < VC2
                    weight4(i,j) = hungry(i) * nPop / sumHungry * rand();
                else
                    weight4(i,j) = 1;
                end
            end
        end

        % Update positions
        shrink = 2 * (1 - it / maxItr);
        for i = 1:nPop
            if rand < VC2
                X(i,:) = X(i,:) .* (1 + randn(1, dim));
            else
                A = randi([1, count]);
                for j = 1:dim
                    r = rand();
                    vb = 2 * shrink * r - shrink;
                    if r > VC1(i)
                        X(i,j) = weight4(i,j)*tempPosition(A,j) + vb*weight3(i,j)*abs(tempPosition(A,j) - X(i,j));
                    else
                        X(i,j) = weight4(i,j)*tempPosition(A,j) - vb*weight3(i,j)*abs(tempPosition(A,j) - X(i,j));
                    end
                end
            end
        end

        convergenceCurve(it) = Destination_fitness;
    end

    bestFitness = Destination_fitness;

end

%% Helper Function: Population Initialization
function Positions = initializePopulation(N, dim, ub, lb)
    Positions = zeros(N, dim);
    if numel(ub) == 1
        Positions = rand(N, dim) * (ub - lb) + lb;
    else
        for i = 1:dim
            Positions(:,i) = rand(N,1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end
