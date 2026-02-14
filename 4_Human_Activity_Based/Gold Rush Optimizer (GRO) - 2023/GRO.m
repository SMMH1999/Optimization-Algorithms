function [bestFitness, bestPosition, convergenceCurve] = GRO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Gold Rush Optimizer (GRO)
    % =========================================================================
    % Author: Kamran Zolfi
    % Reference:
    %   Zolfi, K. (2023). Gold Rush Optimizer: A new population-based
    %   metaheuristic algorithm. Operations Research and Decisions, 33(1).
    %   DOI: 10.37190/ord230108
    %
    % Description:
    %   Gold Rush Optimizer (GRO) is a population-based metaheuristic algorithm
    %   inspired by gold prospectors’ behaviors including migration, mining,
    %   and collaboration strategies.
    %
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Problem dimension (integer)
    %   nPop      - Number of search agents (population size)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found (scalar)
    %   bestPosition     - Best solution found (1×dim vector)
    %   convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % Tunable Parameters:
    %   sigma_initial - Initial exploration coefficient (default = 2)
    %   sigma_final   - Final exploitation coefficient (default = 1/maxItr)
    %
    % =========================================================================
    % Initialization
    % =========================================================================

    % Ensure row vector bounds
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    sigma_initial = 2;
    sigma_final   = 1 / maxItr;

    % Initialize population
    Positions = rand(nPop, dim) .* (ub - lb) + lb;
    X_NEW     = Positions;

    Fit       = inf(1, nPop);
    Fit_NEW   = inf(1, nPop);

    bestFitness  = inf;
    bestPosition = zeros(1, dim);

    convergenceCurve = zeros(1, maxItr);

    % =========================================================================
    % Initial Fitness Evaluation
    % =========================================================================
    for i = 1:nPop
        Fit(i) = objFun(Positions(i,:));
        if Fit(i) < bestFitness
            bestFitness  = Fit(i);
            bestPosition = Positions(i,:);
        end
    end

    % =========================================================================
    % Main Loop
    % =========================================================================
    for t = 1:maxItr

        % Update l1 and l2 (control parameters)
        if maxItr > 1
            l2 = ((maxItr - t)/(maxItr - 1))^2 * (sigma_initial - sigma_final) + sigma_final;
            l1 = ((maxItr - t)/(maxItr - 1))^1 * (sigma_initial - sigma_final) + sigma_final;
        else
            l1 = sigma_final;
            l2 = sigma_final;
        end

        % Evaluate new positions and update
        for i = 1:nPop

            Fit_NEW(i) = objFun(X_NEW(i,:));

            if Fit_NEW(i) < Fit(i)
                Fit(i) = Fit_NEW(i);
                Positions(i,:) = X_NEW(i,:);
            end

            if Fit(i) < bestFitness
                bestFitness  = Fit(i);
                bestPosition = Positions(i,:);
            end
        end

        % Generate next positions
        for i = 1:nPop

            idx = 1:nPop;
            idx(i) = [];
            coworkers = idx(randperm(nPop-1,2));
            digger1 = coworkers(1);
            digger2 = coworkers(2);

            m = rand;

            if m < 1/3
                % Collaboration
                r1 = rand(1,dim);
                D3 = Positions(digger2,:) - Positions(digger1,:);
                X_NEW(i,:) = Positions(i,:) + r1 .* D3;

            elseif m < 2/3
                % Mining
                r1 = rand(1,dim);
                A2 = 2*l2*r1 - l2;
                D2 = Positions(i,:) - Positions(digger1,:);
                X_NEW(i,:) = Positions(digger1,:) + A2 .* D2;

            else
                % Migration
                r1 = rand(1,dim);
                r2 = rand(1,dim);
                C1 = 2 * r2;
                A1 = 1 + l1 * (r1 - 0.5);
                D1 = C1 .* bestPosition - Positions(i,:);
                X_NEW(i,:) = Positions(i,:) + A1 .* D1;
            end

            % Boundary Control
            X_NEW(i,:) = max(X_NEW(i,:), lb);
            X_NEW(i,:) = min(X_NEW(i,:), ub);

        end

        % Store convergence
        convergenceCurve(t) = bestFitness;

    end

end
