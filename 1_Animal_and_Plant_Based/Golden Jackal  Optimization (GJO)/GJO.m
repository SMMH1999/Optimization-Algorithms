function [bestFitness, bestPosition, convergenceCurve] = GJO(lb, ub, dim, nPop, maxItr, objFun)
    % GOLDEN JACKAL OPTIMIZATION (GJO)
    %__________________________________________________________________________
    % Developed by: Nitish Chopra & Muhammad Mohsin Ansari
    % Main paper: Chopra, N., Ansari, M.M., "Golden Jackal Optimization: A
    % Novel Nature-Inspired Optimizer for Engineering Applications," Expert
    % Systems with Applications, 2022.
    % DOI: https://doi.org/10.1016/j.eswa.2022.116924
    %__________________________________________________________________________
    %
    % INPUTS:
    %   lb         - Lower bound (scalar or 1xdim vector)
    %   ub         - Upper bound (scalar or 1xdim vector)
    %   dim        - Problem dimensionality
    %   nPop       - Number of search agents
    %   maxItr     - Maximum number of iterations
    %   objFun     - Handle to objective function, f(x)
    %
    % OUTPUTS:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Position of the best fitness
    %   convergenceCurve - Best fitness at each iteration
    %
    % TUNABLE PARAMETERS:
    %   None beyond input parameters; algorithm automatically balances
    %   exploration/exploitation using E1 (evading energy) and Levy flights.

    %% Initialization
    bestPosition = zeros(1, dim);        % Male jackal
    bestFitness = inf;

    Female_pos = zeros(1, dim);          % Female jackal
    Female_score = inf;

    % Initialize population
    Positions = initializePopulation(nPop, dim, ub, lb);
    convergenceCurve = zeros(1, maxItr);

    %% Main loop
    for t = 1:maxItr
        % Evaluate fitness and update Male & Female Jackals
        for i = 1:nPop
            % Boundary check
            Positions(i,:) = max(min(Positions(i,:), ub), lb);

            fitness = objFun(Positions(i,:));

            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = Positions(i,:);
            elseif fitness > bestFitness && fitness < Female_score
                Female_score = fitness;
                Female_pos = Positions(i,:);
            end
        end

        % Update exploration/exploitation parameters
        E1 = 1.5 * (1 - t/maxItr);
        RL = 0.05 * levyFlight(nPop, dim, 1.5);

        % Update positions
        for i = 1:nPop
            for j = 1:dim
                r1 = rand();
                E0 = 2*r1 - 1;
                E = E1 * E0; % Evading energy

                if abs(E) < 1
                    % Exploitation
                    D_m = abs(RL(i,j)*bestPosition(j) - Positions(i,j));
                    Male_pos = bestPosition(j) - E * D_m;
                    D_f = abs(RL(i,j)*Female_pos(j) - Positions(i,j));
                    Female_new = Female_pos(j) - E * D_f;
                else
                    % Exploration
                    D_m = abs(bestPosition(j) - RL(i,j)*Positions(i,j));
                    Male_pos = bestPosition(j) - E * D_m;
                    D_f = abs(Female_pos(j) - RL(i,j)*Positions(i,j));
                    Female_new = Female_pos(j) - E * D_f;
                end
                % Update position
                Positions(i,j) = (Male_pos + Female_new)/2;
            end
        end

        % Store convergence
        convergenceCurve(t) = bestFitness;
    end

end

%% ---------------------- Helper Functions ---------------------------

function X = initializePopulation(nAgents, dim, ub, lb)
    % Initializes population within bounds
    Boundary_no = numel(ub);
    X = zeros(nAgents, dim);

    if Boundary_no == 1
        X = rand(nAgents, dim) * (ub - lb) + lb;
    else
        for i = 1:dim
            X(:,i) = rand(nAgents,1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end

function z = levyFlight(n, m, beta)
    % Generates Levy flight steps for exploration
    num = gamma(1+beta) * sin(pi*beta/2);
    den = gamma((1+beta)/2) * beta * 2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta);

    u = randn(n, m) * sigma_u;
    v = randn(n, m);
    z = u ./ abs(v).^(1/beta);
end
