function [bestFitness, bestPosition, convergenceCurve] = HBA(lb, ub, dim, nPop, maxItr, objFun)
    %==========================================================================
    % Honey Badger Algorithm (HBA)
    % Author: Hashim et al., 2021
    %
    % Description:
    %   A bio-inspired metaheuristic optimizer based on the foraging behavior
    %   of honey badgers, capable of solving continuous optimization problems.
    %
    % Inputs:
    %   lb       - Lower bound (scalar or 1 x dim vector)
    %   ub       - Upper bound (scalar or 1 x dim vector)
    %   dim      - Number of decision variables (integer)
    %   nPop     - Population size (integer)
    %   maxItr   - Maximum number of iterations (integer)
    %   objFun   - Handle to objective function (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best solution vector found (1 x dim)
    %   convergenceCurve - Best fitness at each iteration (1 x maxItr)
    %
    % Tunable Parameters:
    %   beta  = 6      -> Food attraction factor
    %   C     = 2      -> Density factor constant
    %==========================================================================

    %% Parameters
    beta = 6;               % Food attraction factor
    C = 2;                  % Density factor constant
    vec_flag = [1, -1];     % Direction flag

    %% Initialization
    X = initializePopulation(nPop, dim, ub, lb);         % Population
    fitness = evaluateFitness(objFun, X);               % Fitness evaluation
    [bestFitness, gbestIdx] = min(fitness);
    bestPosition = X(gbestIdx, :);                      % Global best
    convergenceCurve = zeros(1, maxItr);               % Convergence tracking

    %% Main Loop
    for t = 1:maxItr
        alpha = C * exp(-t / maxItr);                  % Density factor
        I = calculateIntensity(nPop, X, bestPosition); % Intensity vector

        for i = 1:nPop
            r = rand();
            F = vec_flag(randi(2));
            Xnew = zeros(1, dim);

            for j = 1:dim
                di = bestPosition(j) - X(i,j);

                if r < 0.5
                    r3 = rand(); r4 = rand(); r5 = rand();
                    Xnew(j) = bestPosition(j) + F*beta*I(i)*bestPosition(j) + ...
                        F*r3*alpha*di*abs(cos(2*pi*r4)*(1 - cos(2*pi*r5)));
                else
                    r7 = rand();
                    Xnew(j) = bestPosition(j) + F*r7*alpha*di;
                end
            end

            % Boundary control
            Xnew = max(min(Xnew, ub), lb);

            % Update if improved
            tempFitness = objFun(Xnew);
            if tempFitness < fitness(i)
                fitness(i) = tempFitness;
                X(i,:) = Xnew;
            end
        end

        % Population-level boundary enforcement
        X = max(min(X, ub), lb);

        % Update global best
        [currentBest, idx] = min(fitness);
        convergenceCurve(t) = currentBest;
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(idx, :);
        end
    end

end

%%=========================================================================
%% Helper Functions
%%=========================================================================

function X = initializePopulation(N, dim, up, down)
    % Initialize population with scalar or vector bounds
    X = zeros(N, dim);
    if isscalar(up)
        X = rand(N, dim) .* (up - down) + down;
    else
        for j = 1:dim
            X(:, j) = rand(N, 1) * (up(j) - down(j)) + down(j);
        end
    end
end

function fitness = evaluateFitness(objFun, X)
    % Evaluate objective function for entire population
    N = size(X,1);
    fitness = zeros(N,1);
    for i = 1:N
        fitness(i) = objFun(X(i,:));
    end
end

function I = calculateIntensity(N, X, Xprey)
    % Calculate intensity vector based on distance and neighbor influence
    di = zeros(1,N); S = zeros(1,N);
    for i = 1:N-1
        di(i) = norm(X(i,:) - Xprey + eps)^2;
        S(i) = norm(X(i,:) - X(i+1,:) + eps)^2;
    end
    di(N) = norm(X(N,:) - Xprey + eps)^2;
    S(N) = norm(X(N,:) - X(1,:) + eps)^2;

    I = zeros(1,N);
    for i = 1:N
        I(i) = rand() * S(i) / (4 * pi * di(i));
    end
end
