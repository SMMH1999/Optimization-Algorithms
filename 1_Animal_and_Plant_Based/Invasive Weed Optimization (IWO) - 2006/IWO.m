function [bestFitness, bestPosition, convergenceCurve] = IWO(lb, ub, dim, nPop, maxItr, objFun)
    %--------------------------------------------------------------------------
    % Invasive Weed Optimization (IWO)
    %
    % Author: S. Mostapha Kalami Heris (Yarpiz Team)
    % Description: Implementation of the IWO algorithm for minimization problems
    %
    % INPUTS:
    %   lb        - Lower bound (scalar or 1xdim vector)
    %   ub        - Upper bound (scalar or 1xdim vector)
    %   dim       - Number of decision variables
    %   nPop      - Maximum population size
    %   maxItr    - Maximum number of iterations
    %   objFun    - Handle to objective function, e.g., @(x) Sphere(x)
    %
    % OUTPUTS:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Position vector of the best solution
    %   convergenceCurve - Best cost per iteration
    %
    % Tunable Parameters (set inside function):
    %   nPop0      - Initial population size
    %   Smin       - Minimum number of seeds
    %   Smax       - Maximum number of seeds
    %   sigma_initial - Initial standard deviation
    %   sigma_final   - Final standard deviation
    %   Exponent      - Variance reduction exponent
    %--------------------------------------------------------------------------

    %% Parameters
    nPop0 = max(5, floor(nPop/3));  % Initial population size
    Smin = 0;                        % Minimum seeds
    Smax = 5;                        % Maximum seeds
    sigma_initial = 0.5;             % Initial std dev
    sigma_final = 0.001;             % Final std dev
    Exponent = 2;                     % Variance reduction exponent

    %% Ensure bounds are vectors
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end

    VarSize = [1 dim];                % Decision variable size

    %% Empty Individual Template
    empty_plant.Position = [];
    empty_plant.Cost = [];

    %% Initialize Population
    pop = repmat(empty_plant, nPop0, 1);
    for i = 1:nPop0
        pop(i).Position = unifrnd(lb, ub, VarSize);
        pop(i).Cost = objFun(pop(i).Position);
    end

    %% Initialize Best Cost History
    convergenceCurve = zeros(maxItr,1);

    %% Main IWO Loop
    for it = 1:maxItr

        % Update Standard Deviation
        sigma = ((maxItr - it)/(maxItr - 1))^Exponent * (sigma_initial - sigma_final) + sigma_final;

        % Get Costs
        Costs = [pop.Cost];
        BestCost = min(Costs);
        WorstCost = max(Costs);

        % Generate Offsprings
        newpop = [];
        for i = 1:numel(pop)
            if BestCost == WorstCost
                ratio = 0;
            else
                ratio = (pop(i).Cost - WorstCost)/(BestCost - WorstCost);
            end
            S = floor(Smin + (Smax - Smin)*ratio);

            for j = 1:S
                newsol = empty_plant;
                newsol.Position = pop(i).Position + sigma*randn(VarSize);
                % Apply bounds
                newsol.Position = max(newsol.Position, lb);
                newsol.Position = min(newsol.Position, ub);
                newsol.Cost = objFun(newsol.Position);
                newpop = [newpop; newsol]; %#ok
            end
        end

        % Merge and Sort
        pop = [pop; newpop];
        [~, SortOrder] = sort([pop.Cost]);
        pop = pop(SortOrder);

        % Competitive Exclusion
        if numel(pop) > nPop
            pop = pop(1:nPop);
        end

        % Update Best Solution
        BestSol = pop(1);
        convergenceCurve(it) = BestSol.Cost;

    end

    %% Return Results
    bestFitness = BestSol.Cost;
    bestPosition = BestSol.Position;

end
