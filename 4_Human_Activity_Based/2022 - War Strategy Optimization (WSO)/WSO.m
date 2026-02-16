function [bestFitness, bestPosition, convergenceCurve] = WSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % War Strategy Optimization (WSO)
    % =========================================================================
    % Author: Dr. Tummala S. L. V. Ayyarao
    % Reference:
    % Ayyarao, T.S.L.V., et al., "War Strategy Optimization Algorithm:
    % A New Effective Metaheuristic Algorithm for Global Optimization,"
    % IEEE Access, 2022.
    %
    % Description:
    % War Strategy Optimization (WSO) is a population-based metaheuristic
    % inspired by strategic war mechanisms including leadership (King),
    % cooperation, and adaptive attack strategies.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of soldiers)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Best solution vector found
    %   convergenceCurve - Best fitness at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   R  - Probability control parameter (default = 0.1)
    %
    % Notes:
    %   • Minimization problem assumed.
    %   • Boundary constraints are strictly enforced.
    %   • No display or external output is generated.
    % =========================================================================

    %% --------------------------- Initialization -----------------------------

    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    Positions = rand(nPop, dim) .* (ub - lb) + lb;

    bestPosition = zeros(1, dim);
    bestFitness = inf;

    fitness = zeros(1, nPop);
    for i = 1:nPop
        fitness(i) = objFun(Positions(i,:));
        if fitness(i) < bestFitness
            bestFitness = fitness(i);
            bestPosition = Positions(i,:);
        end
    end

    fitness_old = fitness;
    fitness_new = inf(1, nPop);

    Positions_new = zeros(size(Positions));

    W1 = 2 * ones(1, nPop);
    Wg = zeros(1, nPop);

    convergenceCurve = zeros(1, maxItr);
    convergenceCurve(1) = bestFitness;

    R = 0.1;

    %% ------------------------------ Main Loop -------------------------------

    for t = 2:maxItr

        [~, sortedIdx] = sort(fitness_old);
        Co = Positions(sortedIdx(2), :);  % Second best

        permIdx = randperm(nPop);

        for i = 1:nPop

            RR = rand;

            if RR < R
                D_V = 2 * RR * (bestPosition - Positions(permIdx(i),:)) + ...
                    W1(i) * rand * (Co - Positions(i,:));
            else
                D_V = 2 * RR * (Co - bestPosition) + ...
                    rand * (W1(i) * bestPosition - Positions(i,:));
            end

            Positions_new(i,:) = Positions(i,:) + D_V;

            % Boundary control
            Positions_new(i,:) = max(Positions_new(i,:), lb);
            Positions_new(i,:) = min(Positions_new(i,:), ub);

            fitness_new(i) = objFun(Positions_new(i,:));

            % Global best update
            if fitness_new(i) < bestFitness
                bestFitness = fitness_new(i);
                bestPosition = Positions_new(i,:);
            end

            % Greedy selection
            if fitness_new(i) < fitness_old(i)
                Positions(i,:) = Positions_new(i,:);
                fitness_old(i) = fitness_new(i);
                Wg(i) = Wg(i) + 1;
                W1(i) = W1(i) * (1 - Wg(i)/maxItr)^2;
            end
        end

        % Re-initialize worst individual (early exploration phase)
        if t < 1000
            [~, worstIdx] = max(fitness_old);
            Positions(worstIdx,:) = rand(1, dim) .* (ub - lb) + lb;
            fitness_old(worstIdx) = objFun(Positions(worstIdx,:));

            if fitness_old(worstIdx) < bestFitness
                bestFitness = fitness_old(worstIdx);
                bestPosition = Positions(worstIdx,:);
            end
        end

        convergenceCurve(t) = bestFitness;

    end

end
