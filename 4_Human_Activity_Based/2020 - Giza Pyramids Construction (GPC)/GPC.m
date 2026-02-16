function [bestFitness, bestPosition, convergenceCurve] = GPC(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Giza Pyramids Construction (GPC) Algorithm
    % -------------------------------------------------------------------------
    % Author   : Sasan Harifi
    % Paper    : Giza Pyramids Construction: an ancient-inspired metaheuristic
    %            algorithm for optimization
    % DOI      : http://dx.doi.org/10.1007/s12065-020-00451-3
    %
    % Description:
    % GPC is a population-based metaheuristic inspired by the ancient pyramid
    % construction process. Each agent (worker) updates its position based on
    % gravity-driven motion, ramp mechanics, friction coefficient, and a
    % stochastic substitution strategy.
    %
    % Input Parameters:
    %   lb        : Lower bound (scalar or 1×dim vector)
    %   ub        : Upper bound (scalar or 1×dim vector)
    %   dim       : Number of decision variables
    %   nPop      : Population size (number of workers)
    %   maxItr    : Maximum number of iterations
    %   objFun    : Objective function handle (minimization)
    %
    % Output:
    %   bestFitness      : Best objective function value found
    %   bestPosition     : Best decision vector found
    %   convergenceCurve : Best fitness value at each iteration
    %
    % Tunable Parameters (Original Paper Defaults):
    %   G       : Gravity constant (9.8)
    %   Tetha   : Ramp angle in degrees (14)
    %   MuMin   : Minimum friction coefficient (1)
    %   MuMax   : Maximum friction coefficient (10)
    %   pSS     : Substitution probability (0.5)
    %
    % =========================================================================

    %% Parameter Settings
    G      = 9.8;
    Tetha  = 14;
    MuMin  = 1;
    MuMax  = 10;
    pSS    = 0.5;

    %% Bound Handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% Initialization
    pop = zeros(nPop, dim);
    fitness = inf(nPop, 1);

    for i = 1:nPop
        pop(i,:) = lb + rand(1, dim) .* (ub - lb);
        fitness(i) = objFun(pop(i,:));
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = pop(idx,:);

    convergenceCurve = zeros(maxItr,1);

    %% Main Loop
    for t = 1:maxItr

        newPop = zeros(nPop, dim);
        newFitness = inf(nPop,1);

        for i = 1:nPop

            V0 = rand;
            Mu = MuMin + (MuMax - MuMin) * rand;

            d = (V0^2) / ((2*G) * (sind(Tetha) + (Mu * cosd(Tetha))));
            x = (V0^2) / ((2*G) * (sind(Tetha)));

            epsilon = (rand(1,dim)-0.5) .* (ub - lb);

            newPos = (pop(i,:) + d) .* (x .* epsilon);

            % Boundary Control
            newPos = max(newPos, lb);
            newPos = min(newPos, ub);

            % Substitution Strategy
            z = pop(i,:);
            k0 = randi(dim);
            for k = 1:dim
                if k == k0 || rand <= pSS
                    z(k) = newPos(k);
                end
            end

            newFitness_i = objFun(z);

            newPop(i,:) = z;
            newFitness(i) = newFitness_i;

            if newFitness_i < bestFitness
                bestFitness = newFitness_i;
                bestPosition = z;
            end
        end

        % Merge
        pop = [pop; newPop];
        fitness = [fitness; newFitness];

        % Sort
        [fitness, sortIdx] = sort(fitness);
        pop = pop(sortIdx,:);

        % Truncate
        pop = pop(1:nPop,:);
        fitness = fitness(1:nPop);

        % Update Global Best
        if fitness(1) < bestFitness
            bestFitness = fitness(1);
            bestPosition = pop(1,:);
        end

        % Convergence Tracking
        convergenceCurve(t) = bestFitness;

    end

end
