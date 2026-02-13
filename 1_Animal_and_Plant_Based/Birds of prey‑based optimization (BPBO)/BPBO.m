function [bestFitness, bestPosition, convergenceCurve] = BPBO(lb, ub, dim, nPop, maxItr, objFun)
    %% Birds of Prey-Based Optimization (BPBO)
    % Author: Mojtaba Ghasemi
    % Co-Author: Nima Khodadadi, University of California Berkeley
    % Email: Nimakhan@berkeley.edu
    % Homepage: https://nimakhodadadi.com
    %
    % Description:
    %   BPBO is a metaheuristic optimization algorithm inspired by birds of prey.
    %   This function minimizes the given objective function 'objFun'.
    %
    % Inputs:
    %   lb        : [1 x dim] Lower bounds of variables (scalar or vector)
    %   ub        : [1 x dim] Upper bounds of variables (scalar or vector)
    %   dim       : Number of decision variables
    %   nPop      : Population size
    %   maxItr    : Maximum number of iterations
    %   objFun    : Handle to the objective function, e.g., @(x) myFunc(x)
    %
    % Outputs:
    %   bestFitness       : Best objective function value found
    %   bestPosition      : Decision variable vector corresponding to bestFitness
    %   convergenceCurve  : Vector of best fitness at each iteration
    %
    % Tunable Parameters:
    %   Pi = 0.7           : Probability factor controlling exploration vs exploitation

    %% Parameter
    Pi = 0.7;

    %% Initialization

    % Empty individual structure
    empty_individual.Position = [];
    empty_individual.Cost = [];

    % Initialize population
    pop = repmat(empty_individual, nPop, 1);

    % Initialize best solution
    Prey.Cost = inf;

    % Initialize population members
    for i = 1:nPop
        pop(i).Position = unifrnd(lb, ub, [1 dim]);
        pop(i).Cost = objFun(pop(i).Position);

        if pop(i).Cost < Prey.Cost
            Prey = pop(i);
        end
    end

    % Initialize convergence record
    convergenceCurve = zeros(maxItr, 1);

    %% Main BPBO Loop
    for it = 1:maxItr

        % Population mean
        Mean = zeros(1, dim);
        for i = 1:nPop
            Mean = Mean + pop(i).Position;
        end
        Mean = Mean / nPop;

        % Select Prey (best individual)
        Prey = pop(1);
        for i = 2:nPop
            if pop(i).Cost < Prey.Cost
                Prey = pop(i);
            end
        end

        % Update each individual
        for i = 1:nPop
            newsol = empty_individual;

            if rand < Pi
                r1 = rand; r2 = rand; r3 = rand;
                if r1 < r2
                    M00 = round(1 + rand);
                    newsol.Position = pop(i).Position + rand([1 dim]).*(Prey.Position - M00*pop(i).Position);
                elseif r1 < r3
                    M01 = round(1 + rand);
                    newsol.Position = Mean + rand([1 dim]).*(Prey.Position - M01*Mean);
                else
                    M02 = round(1 + rand);
                    newsol.Position = pop(i).Position + rand([1 dim]).*(pop(i).Position - M02*pop(nPop).Position);
                end
            else
                newsol.Position = pop(i).Position + rand*unifrnd(lb, ub, [1 dim]);
            end

            % Enforce bounds
            newsol.Position = max(newsol.Position, lb);
            newsol.Position = min(newsol.Position, ub);

            % Evaluate
            newsol.Cost = objFun(newsol.Position);

            % Greedy selection
            if newsol.Cost < pop(i).Cost
                pop(i) = newsol;
                if pop(i).Cost < Prey.Cost
                    Prey = pop(i);
                end
            end
        end

        % Store best cost for this iteration
        convergenceCurve(it) = Prey.Cost;
    end

    %% Outputs
    bestFitness = Prey.Cost;
    bestPosition = Prey.Position;

end
