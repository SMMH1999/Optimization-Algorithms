function [bestFitness, bestPosition, convergenceCurve] = ISA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Interior Search Algorithm (ISA)
    % -------------------------------------------------------------------------
    % Author: Amir H. Gandomi
    % Email: a.h.gandomi@gmail.com
    %
    % References:
    % Gandomi AH. Interior Search Algorithm (ISA): A Novel Approach for
    % Global Optimization. ISA Transactions, 53(4): 1168–1183, 2014.
    % Gandomi AH, Roke DA. Engineering Optimization using Interior Search
    % Algorithm. IEEE Symposium on Computational Intelligence, 2014.
    %
    % -------------------------------------------------------------------------
    % Description:
    % Interior Search Algorithm (ISA) is a population-based metaheuristic
    % inspired by interior design principles including mirror reflection and
    % composition adjustment. The algorithm balances exploration and
    % exploitation using a time-varying control parameter.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Number of decision variables (dimension)
    % nPop      : Population size
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % -------------------------------------------------------------------------
    % Outputs:
    % bestFitness      : Best objective value found
    % bestPosition     : Best decision vector found (1×dim)
    % convergenceCurve : Best fitness value at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Notes:
    % - Minimization problem.
    % - Bound constraints are enforced using ISA evolutionary boundary scheme.
    % - Supports scalar or vector bounds.
    % =========================================================================

    %% ========================= Initialization ===============================
    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end

    ns = lb + rand(nPop, dim) .* (ub - lb);   % Initial population
    x_new = zeros(nPop, dim);

    fvalue = zeros(1, nPop);
    fvalue_new = zeros(1, nPop);

    for i = 1:nPop
        fvalue(i) = objFun(ns(i,:));
    end

    [bestFitness, bestIdx] = min(fvalue);
    bestPosition = ns(bestIdx,:);

    convergenceCurve = zeros(1, maxItr);

    %% =========================== Main Loop ==================================
    for t = 1:maxItr

        alpha = t / maxItr;
        LX = min(ns);
        UX = max(ns);

        [~, bestIdx] = min(fvalue);
        center = ns(bestIdx,:);

        for i = 1:nPop

            if i == bestIdx
                x_new(i,:) = ns(i,:) + randn(1,dim).*(ub-lb)*0.01;
            else
                if rand < alpha
                    beta = rand;
                    mirror = beta*ns(i,:) + (1-beta)*center;
                    x_new(i,:) = 2*mirror - ns(i,:);
                else
                    x_new(i,:) = LX + (UX - LX).*rand(1,dim);
                end
            end

            % ---------------- Boundary Handling ------------------------------
            ns_tmp = x_new(i,:);

            I = ns_tmp < lb;
            if any(I)
                A = rand;
                ns_tmp(I) = A*lb(I) + (1-A)*center(I);
            end

            J = ns_tmp > ub;
            if any(J)
                B = rand;
                ns_tmp(J) = B*ub(J) + (1-B)*center(J);
            end

            x_new(i,:) = ns_tmp;
            % -----------------------------------------------------------------

            fvalue_new(i) = objFun(x_new(i,:));
        end

        % ---------------- Greedy Selection -----------------------------------
        for i = 1:nPop
            if fvalue_new(i) < fvalue(i)
                ns(i,:) = x_new(i,:);
                fvalue(i) = fvalue_new(i);
            end
        end

        % ---------------- Best Update ----------------------------------------
        [currentBest, bestIdx] = min(fvalue);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = ns(bestIdx,:);
        end

        % ---------------- Convergence Curve ----------------------------------
        convergenceCurve(t) = bestFitness;

    end

end
