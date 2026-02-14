function [bestFitness, bestPosition, convergenceCurve] = HBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Heap-Based Optimizer (HBO)
    % =========================================================================
    % Author: Qamar Askari
    % Source Paper:
    % Askari, Q., Saeed, M., & Younas, I. (2020).
    % Heap-based optimizer inspired by corporate rank hierarchy for global optimization.
    % Expert Systems with Applications.
    %
    % Description:
    % Heap-Based Optimizer (HBO) is a population-based metaheuristic inspired
    % by corporate rank hierarchy using a heap data structure. Individuals are
    % organized in a D-ary heap where better solutions rise toward the root.
    % Search agents update their positions based on parent and colleague
    % relationships within the heap.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    % lb        : Lower bound (1×dim or scalar)
    % ub        : Upper bound (1×dim or scalar)
    % dim       : Number of decision variables
    % nPop      : Population size (number of search agents)
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % Outputs:
    % bestFitness      : Best objective function value found
    % bestPosition     : Corresponding decision variables
    % convergenceCurve : Best fitness value at each iteration (1×maxItr)
    %
    % Tunable Parameters (internal):
    % degree : Heap branching factor (default = 3)
    % cycles : Number of cycles controlling gamma dynamics (default = 4)
    %
    % =========================================================================

    %% --------------------------- Parameter Setup ----------------------------
    degree = 3;
    cycles = 4;

    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% --------------------------- Initialization -----------------------------
    Solutions = zeros(nPop, dim);
    for i = 1:dim
        Solutions(:, i) = rand(nPop,1) .* (ub(i) - lb(i)) + lb(i);
    end

    fitnessHeap = inf(nPop, 2);  % Column 1: fitness, Column 2: index
    bestFitness = inf;
    bestPosition = zeros(1, dim);

    %% ----------------------- Initial Fitness Evaluation ---------------------
    for i = 1:nPop
        fitness = objFun(Solutions(i,:));
        fitnessHeap(i,:) = [fitness, i];

        % Heapify upward
        t = i;
        while t > 1
            parentInd = floor((t + 1) / degree);
            if fitnessHeap(t,1) >= fitnessHeap(parentInd,1)
                break;
            else
                temp = fitnessHeap(t,:);
                fitnessHeap(t,:) = fitnessHeap(parentInd,:);
                fitnessHeap(parentInd,:) = temp;
            end
            t = parentInd;
        end

        if fitness < bestFitness
            bestFitness = fitness;
            bestPosition = Solutions(i,:);
        end
    end

    %% ---------------------- Colleague Limits Generation ---------------------
    colleaguesLimits = zeros(nPop,2);
    for c = nPop:-1:1
        hi = ceil((log10(c * degree - c + 1) / log10(degree))) - 1;
        lowerLim = ((degree * degree^(hi-1) - 1)/(degree-1) + 1);
        upperLim = (degree * degree^hi - 1)/(degree-1);
        colleaguesLimits(c,1) = lowerLim;
        colleaguesLimits(c,2) = upperLim;
    end

    %% --------------------------- Main Loop ----------------------------------
    convergenceCurve = zeros(1, maxItr);
    itPerCycle = maxItr / cycles;
    qtrCycle = itPerCycle / 4;

    for t = 1:maxItr

        gamma = mod(t, itPerCycle) / qtrCycle;
        gamma = abs(2 - gamma);

        for c = nPop:-1:2

            parentInd = floor((c + 1) / degree);

            curIdx = fitnessHeap(c,2);
            parentIdx = fitnessHeap(parentInd,2);

            curSol = Solutions(curIdx,:);
            parentSol = Solutions(parentIdx,:);

            % Select colleague
            if colleaguesLimits(c,2) > nPop
                colleaguesLimits(c,2) = nPop;
            end

            colleagueInd = c;
            while colleagueInd == c
                colleagueInd = randi([colleaguesLimits(c,1), colleaguesLimits(c,2)]);
            end

            colleagueIdx = fitnessHeap(colleagueInd,2);
            colleagueSol = Solutions(colleagueIdx,:);

            % Position update
            for j = 1:dim
                p1 = (1 - t/maxItr);
                p2 = p1 + (1 - p1)/2;
                r = rand();
                rn = 2*rand() - 1;

                if r < p1
                    continue;
                elseif r < p2
                    D = abs(parentSol(j) - curSol(j));
                    curSol(j) = parentSol(j) + rn * gamma * D;
                else
                    if fitnessHeap(colleagueInd,1) < fitnessHeap(c,1)
                        D = abs(colleagueSol(j) - curSol(j));
                        curSol(j) = colleagueSol(j) + rn * gamma * D;
                    else
                        D = abs(colleagueSol(j) - curSol(j));
                        curSol(j) = curSol(j) + rn * gamma * D;
                    end
                end
            end

            %% --------------------- Boundary Control --------------------------
            curSol = max(curSol, lb);
            curSol = min(curSol, ub);

            %% ---------------------- Fitness Evaluation -----------------------
            newFitness = objFun(curSol);

            if newFitness < fitnessHeap(c,1)
                fitnessHeap(c,1) = newFitness;
                Solutions(curIdx,:) = curSol;
            end

            if newFitness < bestFitness
                bestFitness = newFitness;
                bestPosition = curSol;
            end

            %% ------------------------- Heapify -------------------------------
            tempIndex = c;
            while tempIndex > 1
                parentInd = floor((tempIndex + 1) / degree);
                if fitnessHeap(tempIndex,1) >= fitnessHeap(parentInd,1)
                    break;
                else
                    temp = fitnessHeap(tempIndex,:);
                    fitnessHeap(tempIndex,:) = fitnessHeap(parentInd,:);
                    fitnessHeap(parentInd,:) = temp;
                end
                tempIndex = parentInd;
            end
        end

        %% -------------------- Convergence Curve ------------------------------
        convergenceCurve(t) = bestFitness;
    end

end
