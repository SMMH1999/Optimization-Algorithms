function [bestFitness, bestPosition, convergenceCurve] = DOA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Dream Optimization Algorithm (DOA)
    % =========================================================================
    % Author: (Based on original DOA source code)
    %
    % Description:
    % DOA is a population-based metaheuristic optimization algorithm that
    % divides the population into 5 groups during the exploration phase and
    % applies memory, forgetting, and supplementation strategies. The first
    % 90% of iterations focus on exploration, while the remaining 10% perform
    % exploitation around the global best solution.
    %
    % =========================================================================
    % Input Parameters:
    % lb        - Lower bound (scalar or 1×dim vector)
    % ub        - Upper bound (scalar or 1×dim vector)
    % dim       - Number of decision variables (problem dimension)
    % nPop      - Population size (number of search agents)
    % maxItr    - Maximum number of iterations
    % objFun    - Objective function handle (minimization)
    %
    % =========================================================================
    % Output Parameters:
    % bestFitness      - Best fitness value found
    % bestPosition     - Best solution vector found (1×dim)
    % convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % =========================================================================
    % Algorithm Parameters:
    % - Population divided into 5 groups during exploration
    % - Exploration phase: first 90% of iterations
    % - Exploitation phase: last 10% of iterations
    % - Boundary handling depends on problem dimension (>15 or ≤15)
    % =========================================================================

    %% ========================= Initialization ==============================

    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    % Population initialization
    X = rand(nPop, dim) .* (ub - lb) + lb;

    SELECT = 1:nPop;

    bestPosition = zeros(1, dim);
    bestFitness  = inf;

    groupBestPos = zeros(5, dim);
    groupBestFit = inf(5, 1);

    convergenceCurve = zeros(1, maxItr);

    %% =========================== Main Loop =================================

    for t = 1:maxItr

        if t <= floor(0.9 * maxItr)
            %% ======================= Exploration Phase ======================

            for m = 1:5

                k = randi([ceil(dim/8/m), ceil(dim/3/m)]);

                startIdx = floor((m-1)/5 * nPop) + 1;
                endIdx   = floor(m/5 * nPop);

                % Group best update
                for j = startIdx:endIdx
                    fit = objFun(X(j,:));
                    if fit < groupBestFit(m)
                        groupBestFit(m) = fit;
                        groupBestPos(m,:) = X(j,:);
                    end
                end

                % Memory + forgetting/supplementation
                for j = startIdx:endIdx
                    X(j,:) = groupBestPos(m,:);
                    in = randperm(dim, k);

                    if rand < 0.9
                        for h = 1:k
                            idx = in(h);
                            X(j,idx) = X(j,idx) + ...
                                (rand*(ub(idx)-lb(idx)) + lb(idx)) * ...
                                (cos((t + maxItr/10)*pi/maxItr)+1)/2;

                            if X(j,idx) > ub(idx) || X(j,idx) < lb(idx)
                                if dim > 15
                                    select = SELECT;
                                    select(j) = [];
                                    sel = select(randi(nPop-1));
                                    X(j,idx) = X(sel,idx);
                                else
                                    X(j,idx) = rand*(ub(idx)-lb(idx)) + lb(idx);
                                end
                            end
                        end
                    else
                        for h = 1:k
                            idx = in(h);
                            X(j,idx) = X(randi(nPop), idx);
                        end
                    end
                end

                % Global best update
                if groupBestFit(m) < bestFitness
                    bestFitness = groupBestFit(m);
                    bestPosition = groupBestPos(m,:);
                end
            end

        else
            %% ======================= Exploitation Phase =====================

            % Global best update
            for p = 1:nPop
                fit = objFun(X(p,:));
                if fit < bestFitness
                    bestFitness = fit;
                    bestPosition = X(p,:);
                end
            end

            for j = 1:nPop
                km = max(2, ceil(dim/3));
                k = randi([2, km]);

                X(j,:) = bestPosition;
                in = randperm(dim, k);

                for h = 1:k
                    idx = in(h);
                    X(j,idx) = X(j,idx) + ...
                        (rand*(ub(idx)-lb(idx)) + lb(idx)) * ...
                        (cos(t*pi/maxItr)+1)/2;

                    if X(j,idx) > ub(idx) || X(j,idx) < lb(idx)
                        if dim > 15
                            select = SELECT;
                            select(j) = [];
                            sel = select(randi(nPop-1));
                            X(j,idx) = X(sel,idx);
                        else
                            X(j,idx) = rand*(ub(idx)-lb(idx)) + lb(idx);
                        end
                    end
                end
            end
        end

        %% ====================== Convergence Curve ===========================
        convergenceCurve(t) = bestFitness;

    end

end
