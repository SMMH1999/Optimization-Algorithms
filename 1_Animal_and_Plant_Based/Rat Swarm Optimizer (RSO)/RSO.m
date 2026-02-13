function [bestFitness, bestPosition, convergenceCurve] = RSO(lb, ub, dim, nPop, maxItr, objFun)
    %_________________________________________________________________________%
    %  Rat Swarm Optimizer (RSO)                                              %
    %                                                                         %
    %  Developed in MATLAB R2019b                                             %
    %                                                                         %
    %  Designed and Developed: Dr. Gaurav Dhiman                              %
    %  E-Mail: gdhiman0001@gmail.com, gaurav.dhiman@ieee.org                  %
    %  Homepage: http://www.dhimangaurav.com                                  %
    %                                                                         %
    %  Published paper: G. Dhiman et al.,
    %    "A novel algorithm for global optimization: Rat Swarm Optimizer",
    %    Journal of Ambient Intelligence and Humanized Computing, 2020, DOI:
    %    https://doi.org/10.1007/s12652-020-02580-0                           %
    %_________________________________________________________________________%
    %
    % Inputs:
    %   lb          : Lower bound (scalar or 1 x dim vector)
    %   ub          : Upper bound (scalar or 1 x dim vector)
    %   dim         : Number of dimensions (integer)
    %   nPop        : Number of search agents (integer)
    %   maxItr      : Maximum number of iterations (integer)
    %   objFun      : Handle to objective function
    %
    % Outputs:
    %   bestFitness     : Best fitness value found
    %   bestPosition    : Position vector of best solution
    %   convergenceCurve: Vector of best fitness at each iteration
    %
    % Tunable Parameters:
    %   R : Random initial coefficient controlling exploration (default 1-5)

    %% Initialization
    bestPosition      = zeros(1, dim);
    bestFitness       = inf;
    convergenceCurve  = zeros(1, maxItr);

    Positions = initializePopulation(nPop, dim, ub, lb);

    x = 1; y = 5;
    R = floor((y - x) * rand() + x); % Exploration coefficient
    iteration = 0;

    %% Main Loop
    while iteration < maxItr
        iteration = iteration + 1;

        % Evaluate fitness and update global best
        for i = 1:nPop
            % Enforce bounds
            Positions(i,:) = max(min(Positions(i,:), ub), lb);

            fitness = objFun(Positions(i,:));

            if fitness < bestFitness
                bestFitness   = fitness;
                bestPosition  = Positions(i,:);
            end
        end

        % Update coefficient for current iteration
        A = R - iteration * (R / maxItr);

        % Position update
        for i = 1:nPop
            for j = 1:dim
                C = 2 * rand();
                P_vec = A * Positions(i,j) + abs(C * (bestPosition(j) - Positions(i,j)));
                Positions(i,j) = bestPosition(j) - P_vec;
            end
        end

        convergenceCurve(iteration) = bestFitness;
    end

end

%% Local Initialization Function
function Pos = initializePopulation(nAgents, dim, ub, lb)
    if isscalar(ub)
        Pos = rand(nAgents, dim) .* (ub - lb) + lb;
    else
        Pos = zeros(nAgents, dim);
        for d = 1:dim
            Pos(:,d) = rand(nAgents, 1) .* (ub(d) - lb(d)) + lb(d);
        end
    end
end
