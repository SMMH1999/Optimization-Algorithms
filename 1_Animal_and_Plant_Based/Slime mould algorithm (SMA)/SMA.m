function [bestFitness, bestPosition, convergenceCurve] = SMA(lb, ub, dim, nPop, maxItr, objFun)
    %% Slime Mould Algorithm (SMA)
    % SMA: A physics-inspired optimization algorithm for continuous problems
    %
    % Inputs:
    %   lb        - lower bound (scalar or 1×dim vector)
    %   ub        - upper bound (scalar or 1×dim vector)
    %   dim       - problem dimension (integer)
    %   nPop      - population size (integer)
    %   maxItr    - maximum number of iterations (integer)
    %   objFun    - handle to objective function: fitness = objFun(position)
    %
    % Outputs:
    %   bestFitness      - best fitness value found (scalar)
    %   bestPosition     - best solution vector found (1×dim)
    %   convergenceCurve - best fitness at each iteration (1×maxItr)
    %
    % Author: Shimin Li, Ali Asghar Heidari, Huiling Chen
    % Reference: Li et al., Future Generation Computer Systems, 2020

    %% Initialization
    bestPosition = zeros(1, dim);
    bestFitness = inf;  % for minimization problems
    convergenceCurve = zeros(1, maxItr);

    % Ensure lb and ub are vectors
    lb = ones(1, dim) * lb;
    ub = ones(1, dim) * ub;

    % Initialize population
    X = rand(nPop, dim) .* (ub - lb) + lb;
    weight = ones(nPop, dim);
    AllFitness = inf(nPop, 1);

    z = 0.03;  % random re-initialization probability

    %% Main Loop
    for it = 1:maxItr
        % Evaluate fitness and enforce boundaries
        for i = 1:nPop
            X(i,:) = max(min(X(i,:), ub), lb);
            AllFitness(i) = objFun(X(i,:));
        end

        [sortedFitness, sortIdx] = sort(AllFitness);
        worstFitness = sortedFitness(end);
        currentBest = sortedFitness(1);
        S = currentBest - worstFitness + eps;  % avoid zero division

        % Update slime mould weights
        for i = 1:nPop
            for j = 1:dim
                if i <= nPop/2
                    weight(sortIdx(i), j) = 1 + rand() * log10((currentBest - sortedFitness(i))/S + 1);
                else
                    weight(sortIdx(i), j) = 1 - rand() * log10((currentBest - sortedFitness(i))/S + 1);
                end
            end
        end

        % Update global best
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(sortIdx(1), :);
        end

        % Parameters for position update
        a = atanh(-(it/maxItr)+1);
        b = 1 - it/maxItr;

        % Update positions
        for i = 1:nPop
            if rand < z
                X(i,:) = (ub - lb).*rand(1, dim) + lb;
            else
                p = tanh(abs(AllFitness(i) - bestFitness));
                vb = unifrnd(-a, a, 1, dim);
                vc = unifrnd(-b, b, 1, dim);
                for j = 1:dim
                    r = rand();
                    A = randi([1, nPop]);
                    B = randi([1, nPop]);
                    if r < p
                        X(i,j) = bestPosition(j) + vb(j) * (weight(i,j)*X(A,j) - X(B,j));
                    else
                        X(i,j) = vc(j) * X(i,j);
                    end
                end
            end
        end

        % Store convergence
        convergenceCurve(it) = bestFitness;
    end

end
