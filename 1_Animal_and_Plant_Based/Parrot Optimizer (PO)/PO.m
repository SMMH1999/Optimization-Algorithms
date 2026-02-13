function [bestFitness, bestPosition, convergenceCurve] = PO(lb, ub, dim, nPop, maxItr, objFun)
    %=================================================================================
    % Parrot Optimizer (PO)
    %
    % Inputs:
    %   lb       : lower bound (scalar or 1xdim vector)
    %   ub       : upper bound (scalar or 1xdim vector)
    %   dim      : number of decision variables
    %   nPop     : population size (number of parrots)
    %   maxItr   : maximum number of iterations
    %   objFun   : handle to objective function, e.g., @(x) sum(x.^2)
    %
    % Outputs:
    %   bestFitness      : best objective value found
    %   bestPosition     : position corresponding to bestFitness
    %   convergenceCurve : bestFitness at each iteration
    %
    %=================================================================================

    % Ensure bounds are vectors
    if isscalar(ub)
        ub = ub .* ones(1, dim);
        lb = lb .* ones(1, dim);
    end

    %% Initialization
    X = rand(nPop, dim) .* (ub - lb) + lb; % initial positions
    fitness = zeros(nPop,1);

    for i = 1:nPop
        fitness(i) = objFun(X(i,:));
    end

    [fitness, idx] = sort(fitness);
    X = X(idx, :);
    bestFitness = fitness(1);
    bestPosition = X(1, :);

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for t = 1:maxItr
        alpha = rand()/5;
        sita = rand() * pi;
        X_new = X;

        for j = 1:nPop
            St = randi([1,4]);

            switch St
                case 1 % Foraging behavior
                    X_new(j,:) = (X(j,:) - bestPosition) .* Levy(dim) + rand() * mean(X) * (1 - t/maxItr)^(2*t/maxItr);
                case 2 % Staying behavior
                    X_new(j,:) = X(j,:) + bestPosition .* Levy(dim) + randn() * (1 - t/maxItr);
                case 3 % Communicating behavior
                    H = rand();
                    if H < 0.5
                        X_new(j,:) = X(j,:) + alpha * (1 - t/maxItr) * (X(j,:) - mean(X));
                    else
                        X_new(j,:) = X(j,:) + alpha * (1 - t/maxItr) * exp(-j/(rand()*maxItr));
                    end
                case 4 % Fear of strangers behavior
                    X_new(j,:) = X(j,:) + rand() * cos(pi*t/(2*maxItr)) * (bestPosition - X(j,:)) ...
                        - cos(sita) * (t/maxItr)^(2/maxItr) * (X(j,:) - bestPosition);
            end

            % Boundary control
            X_new(j,:) = min(max(X_new(j,:), lb), ub);

            % Update global best if improved
            fNew = objFun(X_new(j,:));
            if fNew < bestFitness
                bestFitness = fNew;
                bestPosition = X_new(j,:);
            end
        end

        % Update positions and fitness
        X = X_new;
        for j = 1:nPop
            fitness(j) = objFun(X(j,:));
        end

        % Sort population
        [fitness, idx] = sort(fitness);
        X = X(idx, :);

        convergenceCurve(t) = bestFitness;
    end
end

%% Levy flight function
function step = Levy(d)
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    step = u ./ abs(v).^(1/beta);
end