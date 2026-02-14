function [bestFitness, bestPosition, convergenceCurve] = TSA(lb, ub, dim, nPop, maxItr, objFun)
    % TSA: Tunicate Swarm Algorithm (TSA)
    % Developed by Dr. Gaurav Dhiman (http://dhimangaurav.com/)
    %
    % Inputs:
    %   lb       : lower bound (scalar or vector) [1 x dim]
    %   ub       : upper bound (scalar or vector) [1 x dim]
    %   dim      : problem dimension
    %   nPop     : number of search agents
    %   maxItr   : maximum number of iterations
    %   objFun   : handle to the objective function
    %
    % Outputs:
    %   bestFitness      : best objective function value found
    %   bestPosition     : position vector corresponding to bestFitness
    %   convergenceCurve : vector of bestFitness over iterations
    %
    % Tunable Parameters:
    %   The TSA algorithm internally uses random factors, and xr in [1,4]
    %   controls exploration/exploitation behavior.

    % --------------------- Initialization --------------------- %
    bestPosition = zeros(1, dim);
    bestFitness = inf;
    convergenceCurve = zeros(1, maxItr);

    % Initialize population
    Positions = initializePopulation(nPop, dim, ub, lb);

    % ---------------------- Main Loop ------------------------- %
    for t = 1:maxItr

        for i = 1:nPop
            % Boundary enforcement
            Positions(i,:) = max(min(Positions(i,:), ub), lb);

            % Fitness evaluation
            fitness = objFun(Positions(i,:));
            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = Positions(i,:);
            end
        end

        % Random xr factor for exploration
        xr = fix(1 + rand() * (4-1));

        % Position update
        for i = 1:nPop
            for j = 1:dim
                A1 = ((rand() + rand()) - (2*rand())) / xr;
                c2 = rand();
                c3 = rand();

                if i == 1
                    d_pos = abs(bestPosition(j) - c2 * Positions(i,j));
                    if c3 >= 0
                        Positions(i,j) = bestPosition(j) + A1 * d_pos;
                    else
                        Positions(i,j) = bestPosition(j) - A1 * d_pos;
                    end
                else
                    d_pos = abs(bestPosition(j) - c2 * Positions(i,j));
                    if c3 >= 0
                        Pos_temp = bestPosition(j) + A1 * d_pos;
                    else
                        Pos_temp = bestPosition(j) - A1 * d_pos;
                    end
                    Positions(i,j) = (Pos_temp + Positions(i-1,j)) / 2;
                end
            end
        end

        % Store convergence
        convergenceCurve(t) = bestFitness;

    end
end

% ------------------- Helper Function ------------------- %
function Pos = initializePopulation(nAgents, dim, ub, lb)
    % Supports scalar or vector bounds
    if isscalar(ub)
        Pos = rand(nAgents, dim) * (ub - lb) + lb;
    else
        Pos = zeros(nAgents, dim);
        for d = 1:dim
            Pos(:,d) = rand(nAgents,1) * (ub(d) - lb(d)) + lb(d);
        end
    end
end
