function [bestFitness, bestPosition, convergenceCurve] = ABC(lb, ub, dim, nPop, maxItr, objFun)
    % Artificial Bee Colony (ABC) Algorithm
    % Benchmark Version - Single File Implementation
    % Minimization Problem

    %% Parameter Handling

    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% ABC Parameters

    nOnlooker = nPop;                 % Number of onlooker bees
    L = round(0.6 * dim * nPop);      % Abandonment limit
    a = 1;                            % Acceleration coefficient

    %% Initialization

    pop = zeros(nPop, dim);
    fitness = zeros(nPop, 1);
    trial = zeros(nPop, 1);

    for i = 1:nPop
        pop(i,:) = lb + rand(1,dim).*(ub - lb);
        fitness(i) = objFun(pop(i,:));
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = pop(idx,:);

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop

    for t = 1:maxItr

        %% ================= EMPLOYED BEES =================
        for i = 1:nPop

            % Select k ≠ i
            k = randi(nPop);
            while k == i
                k = randi(nPop);
            end

            phi = a * (2*rand(1,dim) - 1);
            newPos = pop(i,:) + phi .* (pop(i,:) - pop(k,:));

            % Boundary control
            newPos = max(newPos, lb);
            newPos = min(newPos, ub);

            newFit = objFun(newPos);

            if newFit <= fitness(i)
                pop(i,:) = newPos;
                fitness(i) = newFit;
                trial(i) = 0;
            else
                trial(i) = trial(i) + 1;
            end
        end

        %% ================= PROBABILITY CALCULATION =================
        meanCost = mean(fitness);
        F = exp(-fitness ./ meanCost);
        P = F ./ sum(F);

        %% ================= ONLOOKER BEES =================
        for m = 1:nOnlooker

            % Roulette Wheel Selection
            r = rand;
            C = cumsum(P);
            i = find(r <= C, 1, 'first');

            % Select k ≠ i
            k = randi(nPop);
            while k == i
                k = randi(nPop);
            end

            phi = a * (2*rand(1,dim) - 1);
            newPos = pop(i,:) + phi .* (pop(i,:) - pop(k,:));

            % Boundary control
            newPos = max(newPos, lb);
            newPos = min(newPos, ub);

            newFit = objFun(newPos);

            if newFit <= fitness(i)
                pop(i,:) = newPos;
                fitness(i) = newFit;
                trial(i) = 0;
            else
                trial(i) = trial(i) + 1;
            end
        end

        %% ================= SCOUT BEES =================
        for i = 1:nPop
            if trial(i) >= L
                pop(i,:) = lb + rand(1,dim).*(ub - lb);
                fitness(i) = objFun(pop(i,:));
                trial(i) = 0;
            end
        end

        %% ================= UPDATE GLOBAL BEST =================
        [currentBest, idx] = min(fitness);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = pop(idx,:);
        end

        convergenceCurve(t) = bestFitness;

    end

end
