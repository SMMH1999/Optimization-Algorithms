function [bestFitness, bestPosition, convergenceCurve] = CO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    %  Cheetah Optimizer (CO)
    %  =========================================================================
    %  Developed in MATLAB R2018b
    %
    %  Original Authors:
    %  M. A. Akbari, M. Zare, R. Azizipanah-Abarghooei,
    %  S. Mirjalili, M. Dericheh
    %
    %  Title:
    %  The Cheetah Optimizer: A Nature-Inspired Metaheuristic Algorithm
    %  for Large-Scale Optimization Problems
    %  Scientific Reports
    %
    % -------------------------------------------------------------------------
    %  Description:
    %  CO simulates cooperative hunting behavior of cheetahs including
    %  searching, attacking, sit-and-wait strategies, leader update,
    %  prey escape, and home-return mechanism.
    %
    % -------------------------------------------------------------------------
    %  Inputs:
    %    lb        : Lower bound (scalar or 1×dim vector)
    %    ub        : Upper bound (scalar or 1×dim vector)
    %    dim       : Number of decision variables
    %    nPop      : Population size (number of cheetahs)
    %    maxItr    : Maximum number of function evaluations
    %    objFun    : Objective function handle (minimization)
    %
    % -------------------------------------------------------------------------
    %  Outputs:
    %    bestFitness      : Best fitness value found
    %    bestPosition     : Best solution vector
    %    convergenceCurve : Best fitness per iteration (prey fitness)
    %
    % -------------------------------------------------------------------------
    %  Internal Parameters:
    %    m  : Number of search agents in hunting group (fixed = 2)
    %    T  : Hunting time = ceil(dim/10)*60
    % =========================================================================

    %% ============================ Initialization ============================
    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    end

    m = 2;                                   % Group size
    T = ceil(dim/10) * 60;                   % Hunting time
    MaxFEs = maxItr;                         % Maximum evaluations

    empty_individual.Position = [];
    empty_individual.Cost = [];
    pop = repmat(empty_individual, nPop, 1);

    BestSol.Cost = inf;

    % Initial population
    for i = 1:nPop
        pop(i).Position = lb + rand(1, dim) .* (ub - lb);
        pop(i).Cost = objFun(pop(i).Position);
        if pop(i).Cost < BestSol.Cost
            BestSol = pop(i);
        end
    end

    pop1 = pop;                % Initial home positions
    X_best = BestSol;          % Prey (global best so far)

    BestCost = [];
    Globest = [];

    t = 0;                     % Hunting time counter
    it = 1;                    % Iteration counter
    FEs = nPop;                % Function evaluations counter

    convergenceCurve = zeros(1, MaxFEs);

    %% ============================ Main Loop ================================
    while FEs <= MaxFEs

        i0 = randi(nPop, 1, m);      % Random hunting group

        for k = 1:m

            i = i0(k);

            % Neighbor selection
            if k == m
                a = i0(k-1);
            else
                a = i0(k+1);
            end

            X  = pop(i).Position;
            X1 = pop(a).Position;
            Xb = BestSol.Position;
            Xbest = X_best.Position;

            kk = 0;

            if i <= 2 && t > 2 && t > ceil(0.2*T+1) && ...
                    abs(BestCost(t-2) - BestCost(t-ceil(0.2*T+1))) <= ...
                    0.0001 * Globest(t-1)
                X = X_best.Position;
                kk = 0;
            elseif i == 3
                X = BestSol.Position;
                kk = -0.1 * rand * t / T;
            else
                kk = 0.25;
            end

            Z = X;

            % ================= Strategy Update =================
            r_Hat = randn(1, dim);
            r1    = rand(1, dim);

            if k == 1
                alpha = 0.0001 * t/T .* (ub - lb);
            else
                alpha = 0.0001 * t/T .* abs(Xb - X) + ...
                    0.001 .* round(double(rand(1,dim) > 0.6));
            end

            r = randn(1, dim);
            r_Check = abs(r).^exp(r./2) .* sin(2*pi*r);
            beta = X1 - X;

            h0 = exp(2 - 2*t/T);
            H  = abs(2 .* r1 .* h0 - h0);

            r2 = rand(1, dim);
            r3 = kk + rand(1, dim);
            r4 = 3 * rand(1, dim);

            % Search & Attack
            f1 = find(H > r4);
            Z(f1) = X(f1) + r_Hat(f1).^-1 .* alpha(f1);

            f1 = find(H <= r4);
            Z(f1) = Xbest(f1) + r_Check(f1) .* beta(f1);

            % Sit & Wait
            f1 = find(r2 > r3);
            Z(f1) = X(f1);

            % ================= Boundary Control =================
            idx = Z < lb;
            Z(idx) = lb(idx) + rand(1,sum(idx)) .* (ub(idx)-lb(idx));

            idx = Z > ub;
            Z(idx) = lb(idx) + rand(1,sum(idx)) .* (ub(idx)-lb(idx));

            % ================= Evaluation =================
            NewSol.Position = Z;
            NewSol.Cost = objFun(Z);
            FEs = FEs + 1;

            if NewSol.Cost < pop(i).Cost
                pop(i) = NewSol;
                if pop(i).Cost < BestSol.Cost
                    BestSol = pop(i);
                end
            end

            if FEs > MaxFEs
                break;
            end
        end

        t = t + 1;
        it = it + 1;

        % ================= Leave Prey & Return Home =================
        if t > T && t-round(T)-1 >= 1 && t > 2
            if abs(BestCost(t-1) - BestCost(t-round(T)-1)) <= ...
                    abs(0.01 * BestCost(t-1))

                best = X_best.Position;
                j0 = randi(dim,1,ceil(dim/10*rand));
                best(j0) = lb(j0) + rand(1,length(j0)) .* ...
                    (ub(j0) - lb(j0));

                BestSol.Position = best;
                BestSol.Cost = objFun(best);
                FEs = FEs + 1;

                i0 = randi(nPop,1,nPop);
                pop(i0(nPop-m+1:nPop)) = pop1(i0(1:m));
                pop(i) = X_best;

                t = 1;
            end
        end

        % ================= Global Best Update =================
        if BestSol.Cost < X_best.Cost
            X_best = BestSol;
        end

        BestCost(t) = BestSol.Cost;
        Globest(t)  = X_best.Cost;

        convergenceCurve(min(FEs,MaxFEs)) = X_best.Cost;

    end

    %% ============================ Output ====================================
    bestFitness  = X_best.Cost;
    bestPosition = X_best.Position;

    convergenceCurve = convergenceCurve(1:FEs-1);

end
