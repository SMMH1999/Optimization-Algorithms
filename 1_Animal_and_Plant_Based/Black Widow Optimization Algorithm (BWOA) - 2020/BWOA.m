function [bestFitness, bestPosition, convergenceCurve] = BWOA(lb, ub, dim, nPop, maxItr, objFun)
    %_________________________________________________________________________%
    % Black Widow Optimization Algorithm (BWOA)                                %
    % Developed in MATLAB R2018a (9.4)                                         %
    % Author: Dr. Hernan Peraza-Vazquez                                        %
    % Programmer: Hernan Peraza-Vazquez                                         %
    % e-Mails: hperaza@ipn.mx or hperaza@ieee.org                              %
    % Paper: A Novel Bio-Inspired Algorithm Applied to Selective Harmonic       %
    %        Elimination in a Three-Phase Eleven-Level Inverter                %
    %        https://doi.org/10.1155/2020/8856040                               %
    %_________________________________________________________________________%
    % INPUTS:
    %   lb         : scalar or 1xD vector, lower bounds of search space
    %   ub         : scalar or 1xD vector, upper bounds of search space
    %   dim        : integer, number of decision variables
    %   nPop       : integer, number of search agents (population size)
    %   maxItr     : integer, maximum number of iterations
    %   objFun     : function handle, objective function to minimize
    %
    % OUTPUTS:
    %   bestFitness      : best (minimum) fitness value found
    %   bestPosition     : 1xD vector, position corresponding to bestFitness
    %   convergenceCurve : 1x(maxItr+1) vector, best fitness at each iteration
    %_________________________________________________________________________%

    %% Initialization
    Positions = initializePopulation(nPop, dim, ub, lb);
    Fitness = zeros(nPop, 1);

    for i = 1:nPop
        Fitness(i) = objFun(Positions(i, :));
    end

    [bestFitness, minIdx] = min(Fitness);
    bestPosition = Positions(minIdx, :);

    [vMax, ~] = max(Fitness);
    convergenceCurve = zeros(1, maxItr + 1);
    convergenceCurve(1) = bestFitness;

    pheromone = computePheromone(Fitness, bestFitness, vMax);

    %% Main Loop
    for t = 1:maxItr
        beta = -1 + 2 * rand();       % -1 < beta < 1
        m = 0.4 + 0.5 * rand();       % 0.4 < m < 0.9

        for r = 1:nPop
            P = rand();
            r1 = randi([1, nPop]);

            % Spiral or direct search
            if P >= 0.3
                v = bestPosition - cos(2*pi*beta) * Positions(r, :);
            else
                v = bestPosition - m * Positions(r1, :);
            end

            % Pheromone-based adjustment
            if pheromone(r) <= 0.3
                band = 1;
                while band
                    r1 = randi([1, nPop]);
                    r2 = randi([1, nPop]);
                    if r1 ~= r2
                        band = 0;
                    end
                end
                v = bestPosition + (Positions(r1,:) - ((-1)^getBinary()) * Positions(r2,:)) / 2;
            end

            % Boundary check
            v = max(min(v, ub), lb);

            % Evaluate fitness
            Fnew = objFun(v);

            % Update position if improved
            if Fnew <= Fitness(r)
                Positions(r, :) = v;
                Fitness(r) = Fnew;
            end

            % Update global best
            if Fnew <= bestFitness
                bestPosition = v;
                bestFitness = Fnew;
            end
        end

        % Update pheromones
        [vMax, ~] = max(Fitness);
        pheromone = computePheromone(Fitness, bestFitness, vMax);

        % Record convergence
        convergenceCurve(t+1) = bestFitness;
    end

end

%% ------------------ Helper Functions ------------------ %

function val = getBinary()
    % Algorithm-2 in paper (page 6)
    val = rand() >= 0.5;
end

function pheromone = computePheromone(fit, fmin, fmax)
    % Eq. 12 in the paper
    pheromone = (fmax - fit) ./ (fmax - fmin + eps);
end

function Positions = initializePopulation(nPop, dim, ub, lb)
    % Initialize population within bounds (scalar or vector)
    if isscalar(ub)
        Positions = rand(nPop, dim) .* (ub - lb) + lb;
    else
        Positions = zeros(nPop, dim);
        for i = 1:dim
            Positions(:, i) = rand(nPop, 1) .* (ub(i) - lb(i)) + lb(i);
        end
    end
end
