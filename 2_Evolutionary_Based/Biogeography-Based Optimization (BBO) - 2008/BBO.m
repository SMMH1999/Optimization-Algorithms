function [bestFitness, bestPosition, convergenceCurve] = BBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Biogeography-Based Optimization (BBO)
    % -------------------------------------------------------------------------
    % Algorithm: Biogeography-Based Optimization (BBO)
    % Abbreviation: BBO
    %
    % Developer (Original Version): S. Mostapha Kalami Heris (Yarpiz Team)
    % Refactored: Unified benchmark-style implementation
    %
    % Description:
    % BBO is a population-based evolutionary optimization algorithm inspired
    % by the migration behavior of species between habitats. Solutions are
    % represented as habitats, and information sharing occurs through
    % immigration and emigration operators. Mutation introduces additional
    % diversity.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size (number of habitats)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (maxItr×1)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters (Fixed as in original code):
    %   KeepRate  - Fraction of best habitats preserved each iteration (0.2)
    %   alpha     - Migration step size factor (0.9)
    %   pMutation - Mutation probability per variable (0.1)
    %   sigma     - Mutation standard deviation (2% of search range)
    %
    % =========================================================================

    %% Parameter Handling

    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    lb = lb(:)';
    ub = ub(:)';

    KeepRate  = 0.2;
    nKeep     = round(KeepRate * nPop);
    nNew      = nPop - nKeep;

    mu        = linspace(1, 0, nPop);   % Emigration rates
    lambda    = 1 - mu;                 % Immigration rates

    alpha     = 0.9;
    pMutation = 0.1;
    sigma     = 0.02 * (ub - lb);

    %% Initialization

    popPos = zeros(nPop, dim);
    popFit = zeros(nPop, 1);

    for i = 1:nPop
        popPos(i, :) = lb + rand(1, dim) .* (ub - lb);
        popFit(i)    = objFun(popPos(i, :));
    end

    [popFit, sortIdx] = sort(popFit);
    popPos = popPos(sortIdx, :);

    bestPosition = popPos(1, :);
    bestFitness  = popFit(1);

    convergenceCurve = zeros(maxItr, 1);

    %% Main Loop

    for t = 1:maxItr

        newPopPos = popPos;
        newPopFit = zeros(nPop, 1);

        for i = 1:nPop

            for k = 1:dim

                % Migration
                if rand <= lambda(i)

                    EP = mu;
                    EP(i) = 0;
                    EP = EP / sum(EP);

                    j = rouletteWheelSelection(EP);

                    newPopPos(i, k) = popPos(i, k) + ...
                        alpha * (popPos(j, k) - popPos(i, k));
                end

                % Mutation
                if rand <= pMutation
                    newPopPos(i, k) = newPopPos(i, k) + sigma(k) * randn;
                end

            end

            % Boundary Control
            newPopPos(i, :) = max(newPopPos(i, :), lb);
            newPopPos(i, :) = min(newPopPos(i, :), ub);

            % Fitness Evaluation
            newPopFit(i) = objFun(newPopPos(i, :));
        end

        % Sort New Population
        [newPopFit, sortIdx] = sort(newPopFit);
        newPopPos = newPopPos(sortIdx, :);

        % Create Next Generation
        popPos = [popPos(1:nKeep, :);
            newPopPos(1:nNew, :)];

        popFit = [popFit(1:nKeep);
            newPopFit(1:nNew)];

        % Sort Combined Population
        [popFit, sortIdx] = sort(popFit);
        popPos = popPos(sortIdx, :);

        % Update Global Best
        if popFit(1) < bestFitness
            bestFitness  = popFit(1);
            bestPosition = popPos(1, :);
        end

        convergenceCurve(t) = bestFitness;

    end

end

%% ------------------------------------------------------------------------
% Roulette Wheel Selection (Local Function)
% -------------------------------------------------------------------------
function index = rouletteWheelSelection(prob)

    r = rand;
    c = cumsum(prob);
    index = find(r <= c, 1, 'first');

    if isempty(index)
        index = numel(prob);
    end

end
