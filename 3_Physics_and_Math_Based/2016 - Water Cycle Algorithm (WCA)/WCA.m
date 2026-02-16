function [bestFitness, bestPosition, convergenceCurve] = WCA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Water Cycle Algorithm (WCA)
    % =========================================================================
    % Standard Version of Water Cycle Algorithm (WCA)
    % Reference:
    % H. Eskandar et al., “Water cycle algorithm - a novel metaheuristic
    % optimization method for solving constrained engineering optimization
    % problems”, Computers & Structures, 110-111 (2012) 151-166.
    %
    % This implementation is adapted to unconstrained single-objective
    % minimization problems for benchmark frameworks.
    %
    % INPUTS:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Number of decision variables
    % nPop      : Population size
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % OUTPUTS:
    % bestFitness      : Best fitness value found
    % bestPosition     : Best decision vector found
    % convergenceCurve : Best fitness value at each iteration
    %
    % Tunable Internal Parameters:
    % Nsr   : Number of rivers + sea (default = 4)
    % dmax  : Evaporation condition constant (default = 1e-5)
    %
    % =========================================================================
    % Author: Original WCA by Eskandar et al.
    % Refactored for benchmark use.
    % =========================================================================

    %% --------------------------- Parameters ---------------------------------
    Nsr  = 4;          % Number of rivers + sea
    dmax = 1e-5;       % Evaporation threshold

    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    %% --------------------------- Initialization ------------------------------
    N_stream = nPop - Nsr;

    empty_ind.position = [];
    empty_ind.cost     = [];

    pop = repmat(empty_ind, nPop, 1);

    for i = 1:nPop
        pop(i).position = lb + rand(1, dim) .* (ub - lb);
        pop(i).cost     = objFun(pop(i).position);
    end

    % Sort population (ascending for minimization)
    [~, idx] = sort([pop.cost]);
    pop = pop(idx);

    sea    = pop(1);
    river  = pop(2:Nsr);
    stream = pop(Nsr+1:end);

    %% ---------------- Stream Allocation --------------------------------------
    costs = [sea.cost; [river.cost]'; [stream.cost]'];
    CN = costs - max(costs);
    if sum(abs(CN)) == 0
        NS = round((1/(Nsr))*N_stream) * ones(1, Nsr);
    else
        NS = round(abs(CN(1:Nsr)) / sum(abs(CN(1:Nsr))) * N_stream);
    end

    % Adjust allocation
    while sum(NS) > N_stream
        [~, id] = max(NS);
        NS(id) = NS(id) - 1;
    end
    while sum(NS) < N_stream
        [~, id] = min(NS);
        NS(id) = NS(id) + 1;
    end

    NB = NS(2:end);

    %% ------------------------- Main Loop ------------------------------------
    convergenceCurve = zeros(maxItr,1);

    for t = 1:maxItr

        % -------- Streams move toward Sea --------
        for j = 1:NS(1)
            stream(j).position = stream(j).position + ...
                2*rand(1,dim).*(sea.position - stream(j).position);

            stream(j).position = min(max(stream(j).position, lb), ub);
            stream(j).cost = objFun(stream(j).position);

            if stream(j).cost < sea.cost
                temp = sea;
                sea = stream(j);
                stream(j) = temp;
            end
        end

        % -------- Streams move toward Rivers --------
        for k = 1:Nsr-1
            for j = 1:NB(k)
                idx_stream = j + sum(NS(1:k));
                stream(idx_stream).position = ...
                    stream(idx_stream).position + ...
                    2*rand(1,dim).*(river(k).position - ...
                    stream(idx_stream).position);

                stream(idx_stream).position = ...
                    min(max(stream(idx_stream).position, lb), ub);

                stream(idx_stream).cost = objFun(stream(idx_stream).position);

                if stream(idx_stream).cost < river(k).cost
                    temp = river(k);
                    river(k) = stream(idx_stream);
                    stream(idx_stream) = temp;

                    if river(k).cost < sea.cost
                        temp2 = sea;
                        sea = river(k);
                        river(k) = temp2;
                    end
                end
            end
        end

        % -------- Rivers move toward Sea --------
        for k = 1:Nsr-1
            river(k).position = river(k).position + ...
                2*rand(1,dim).*(sea.position - river(k).position);

            river(k).position = min(max(river(k).position, lb), ub);
            river(k).cost = objFun(river(k).position);

            if river(k).cost < sea.cost
                temp = sea;
                sea = river(k);
                river(k) = temp;
            end
        end

        % -------- Evaporation and Rain --------
        for k = 1:Nsr-1
            if norm(river(k).position - sea.position) < dmax || rand < 0.1
                for j = 1:NB(k)
                    idx_stream = j + sum(NS(1:k));
                    stream(idx_stream).position = ...
                        lb + rand(1,dim).*(ub - lb);
                    stream(idx_stream).cost = ...
                        objFun(stream(idx_stream).position);
                end
            end
        end

        for j = 1:NS(1)
            if norm(stream(j).position - sea.position) < dmax
                stream(j).position = sea.position + ...
                    sqrt(0.1).*randn(1,dim);
                stream(j).position = ...
                    min(max(stream(j).position, lb), ub);
                stream(j).cost = objFun(stream(j).position);
            end
        end

        dmax = dmax - (dmax/maxItr);

        convergenceCurve(t) = sea.cost;
    end

    %% ------------------------- Final Output ----------------------------------
    bestPosition = sea.position;
    bestFitness  = sea.cost;

end
