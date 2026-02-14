function [bestFitness, bestPosition, convergenceCurve] = CBSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Algorithm: Connected Banking Systems Optimizer(CBSO)
    % =========================================================================
    % Author: Unknown (Refactored to Benchmark Format)
    %
    % Description:
    % This algorithm is a population-based metaheuristic optimizer that uses
    % Gaussian random walks, Levy flights, elitism-based greedy selection,
    % and adaptive coefficient control to balance exploration and exploitation.
    % The algorithm is designed for continuous optimization problems.
    %
    % Inputs:
    %   lb        : Lower bound (scalar or 1×dim vector)
    %   ub        : Upper bound (scalar or 1×dim vector)
    %   dim       : Number of decision variables (dimension)
    %   nPop      : Population size
    %   maxItr    : Maximum number of iterations
    %   objFun    : Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      : Best objective function value found
    %   bestPosition     : Best solution vector found (1×dim)
    %   convergenceCurve : Best fitness value at each iteration (1×maxItr)
    %
    % Tunable Parameters:
    %   beta  : Levy distribution parameter (fixed at 1.5)
    %   CF    : Adaptive control factor (linearly decreasing)
    % =========================================================================

    %% =========================== Initialization =============================
    if isscalar(lb)
        lb = repmat(lb, 1, dim);
    end
    if isscalar(ub)
        ub = repmat(ub, 1, dim);
    end

    lbMat = repmat(lb, nPop, 1);
    ubMat = repmat(ub, nPop, 1);

    X = rand(nPop, dim) .* (ubMat - lbMat) + lbMat;

    fitness = inf(nPop,1);
    bestFitness = inf;
    bestPosition = zeros(1,dim);
    convergenceCurve = zeros(1,maxItr);

    %% ========================= Initial Evaluation ===========================
    for i = 1:nPop
        fitness(i) = objFun(X(i,:));
        if fitness(i) < bestFitness
            bestFitness = fitness(i);
            bestPosition = X(i,:);
        end
    end

    fitness_old = fitness;
    X_old = X;

    %% ============================= Main Loop ================================
    for t = 1:maxItr

        E = repmat(bestPosition, nPop, 1);
        CF = rand * (1 - t/maxItr);

        LD = levyFlight(nPop, dim, 1.5);
        ND = randn(nPop, dim);

        for i = 1:nPop
            SEL1 = randi([1 i]);
            SEL2 = randi([1 i]);

            for j = 1:dim
                if t < 0.2 * maxItr

                    X(i,j) = X(i,j) + ND(i,j) * (rand * E(i,j) - X(i,j));

                elseif t <= 0.4 * maxItr

                    if i > nPop/2
                        X(i,j) = E(i,j) + CF * (ND(i,j) * (rand * E(i,j) - X(SEL1,j)));
                    else
                        X(i,j) = X(i,j) + LD(i,j) * (E(i,j) - LD(i,j) * X(SEL2,j));
                    end

                else

                    X(i,j) = E(i,j) + CF * (LD(i,j) * (LD(i,j) * X(SEL1,j) - X(SEL2,j)));

                end
            end
        end

        %% ======================= Boundary Control ===========================
        X = max(X, lbMat);
        X = min(X, ubMat);

        %% ========================= Fitness Update ============================
        for i = 1:nPop
            fitness(i) = objFun(X(i,:));
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = X(i,:);
            end
        end

        %% ===================== Greedy Selection ==============================
        improveMask = fitness_old < fitness;
        X(improveMask,:) = X_old(improveMask,:);
        fitness(improveMask) = fitness_old(improveMask);

        fitness_old = fitness;
        X_old = X;

        %% =================== Additional Perturbation =========================
        if rand < 0.2
            X = X + rand * (CF * lbMat + rand(nPop,dim) .* (ubMat - lbMat));
        else
            r1 = randperm(nPop);
            r2 = randperm(nPop);
            X = X + rand * rand * (X(r1,:) - X(r2,:));
        end

        %% ===================== Store Convergence =============================
        convergenceCurve(t) = bestFitness;

    end

end

%% =========================== Levy Function ===============================
function z = levyFlight(n, m, beta)

    num = gamma(1+beta) * sin(pi*beta/2);
    den = gamma((1+beta)/2) * beta * 2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta);

    u = sigma_u .* randn(n,m);
    v = randn(n,m);

    z = u ./ (abs(v).^(1/beta));

end
