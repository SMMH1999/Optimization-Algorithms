function [bestFitness, bestPosition, convergenceCurve] = CS(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Cuckoo Search (CS)
    % =========================================================================
    % Algorithm Name: Cuckoo Search (CS)
    % Abbreviation  : CS
    %
    % Original Authors:
    %   Xin-She Yang and Suash Deb
    %   University of Cambridge (2008–2009)
    %
    % Description:
    %   Cuckoo Search is a nature-inspired metaheuristic optimization algorithm
    %   based on brood parasitism of some cuckoo species combined with Lévy
    %   flights. Each cuckoo lays one egg at a time and dumps it in a randomly
    %   chosen nest. A fraction (pa) of worst nests are abandoned and replaced.
    %
    %   This implementation preserves the original logic and behavior of the
    %   standard CS algorithm while reformatted into a unified benchmark-style
    %   function.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb      : Lower bound (scalar or 1×dim vector)
    %   ub      : Upper bound (scalar or 1×dim vector)
    %   dim     : Problem dimension (integer)
    %   nPop    : Number of nests (population size)
    %   maxItr  : Maximum number of iterations
    %   objFun  : Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      : Best objective value found
    %   bestPosition     : Best solution vector found (1×dim)
    %   convergenceCurve : Best fitness value at each iteration (maxItr×1)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   pa   : Discovery rate of alien eggs (default = 0.25)
    %   beta : Lévy distribution exponent (default = 3/2)
    %
    % =========================================================================

    %% --------------------------- Initialization -----------------------------

    pa = 0.25;               % Discovery probability
    beta = 3/2;              % Lévy exponent

    % Bound handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    % Initialize nests
    nest = lb + (ub - lb) .* rand(nPop, dim);

    % Evaluate initial fitness
    fitness = zeros(nPop, 1);
    for i = 1:nPop
        fitness(i) = objFun(nest(i, :));
    end

    % Determine initial best
    [bestFitness, idx] = min(fitness);
    bestPosition = nest(idx, :);

    convergenceCurve = zeros(maxItr, 1);

    %% ---------------------------- Main Loop ---------------------------------

    for t = 1:maxItr

        % --------- Generate New Solutions via Lévy Flights -------------------
        newNest = getCuckoos(nest, bestPosition, lb, ub, beta);

        % Evaluate and Greedy Selection
        for i = 1:nPop
            fnew = objFun(newNest(i, :));
            if fnew <= fitness(i)
                fitness(i) = fnew;
                nest(i, :) = newNest(i, :);
            end
        end

        % --------- Discovery and Randomization -------------------------------
        newNest = emptyNests(nest, lb, ub, pa);

        % Evaluate and Greedy Selection
        for i = 1:nPop
            fnew = objFun(newNest(i, :));
            if fnew <= fitness(i)
                fitness(i) = fnew;
                nest(i, :) = newNest(i, :);
            end
        end

        % --------- Update Global Best ----------------------------------------
        [currentBest, idx] = min(fitness);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = nest(idx, :);
        end

        % Store convergence
        convergenceCurve(t) = bestFitness;
    end

end

%% ========================== Local Functions =============================

function nest = getCuckoos(nest, best, lb, ub, beta)

    nPop = size(nest, 1);
    dim  = size(nest, 2);

    sigma = (gamma(1+beta)*sin(pi*beta/2) / ...
        (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

    for j = 1:nPop
        s = nest(j, :);

        u = randn(1, dim) * sigma;
        v = randn(1, dim);
        step = u ./ abs(v).^(1/beta);

        stepsize = 0.01 * step .* (s - best);
        s = s + stepsize .* randn(1, dim);

        % Apply bounds
        s = simpleBounds(s, lb, ub);
        nest(j, :) = s;
    end

end

% -------------------------------------------------------------------------

function newNest = emptyNests(nest, lb, ub, pa)

    [nPop, dim] = size(nest);

    K = rand(nPop, dim) > pa;

    perm1 = randperm(nPop);
    perm2 = randperm(nPop);

    stepsize = rand * (nest(perm1, :) - nest(perm2, :));
    newNest = nest + stepsize .* K;

    for j = 1:nPop
        newNest(j, :) = simpleBounds(newNest(j, :), lb, ub);
    end

end

% -------------------------------------------------------------------------

function s = simpleBounds(s, lb, ub)

    s = max(s, lb);
    s = min(s, ub);

end
