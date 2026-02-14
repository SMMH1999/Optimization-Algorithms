function [Silverback_Score, Silverback, convergence_curve] = GTO(lb, ub, dim, nPop, maxItr, fobj)
    % GTO - Gorilla Troops Optimizer (Benchmark-ready)
    % INPUTS:
    %   lb       - lower bound (scalar or vector)
    %   ub       - upper bound (scalar or vector)
    %   dim      - number of decision variables
    %   nPop     - population size
    %   maxItr   - maximum number of iterations
    %   fobj     - objective function handle
    % OUTPUTS:
    %   Silverback_Score  - best fitness value found
    %   Silverback        - best solution (position vector)
    %   convergence_curve - vector of best fitness per iteration

    %% Initialize Silverback and population
    Silverback = [];
    Silverback_Score = inf;

    X = initializePopulation(nPop, dim, lb, ub);
    Pop_Fit = zeros(nPop, 1);

    for i = 1:nPop
        Pop_Fit(i) = fobj(X(i,:));
        if Pop_Fit(i) < Silverback_Score
            Silverback_Score = Pop_Fit(i);
            Silverback = X(i,:);
        end
    end

    GX = X;

    % Ensure lb and ub are vectors
    lb = repmat(lb, 1, dim);
    ub = repmat(ub, 1, dim);

    % Parameters
    p = 0.03;
    Beta = 3;
    w = 0.8;

    convergence_curve = zeros(maxItr,1);

    %% Main loop
    for It = 1:maxItr
        a = (cos(2*rand)+1)*(1 - It/maxItr);
        C = a*(2*rand-1);

        %% Exploration
        for i = 1:nPop
            if rand < p
                GX(i,:) = lb + (ub-lb).*rand(1,dim);
            else
                if rand >= 0.5
                    Z = unifrnd(-a, a, 1, dim);
                    H = Z .* X(i,:);
                    GX(i,:) = (rand - a)*X(randi([1,nPop]),:) + C.*H;
                else
                    GX(i,:) = X(i,:) - C.*(C*(X(i,:) - GX(randi([1,nPop]),:)) + rand*(X(i,:) - GX(randi([1,nPop]),:)));
                end
            end
        end

        GX = boundaryCheck(GX, lb, ub);

        % Group formation update
        for i = 1:nPop
            newFit = fobj(GX(i,:));
            if newFit < Pop_Fit(i)
                Pop_Fit(i) = newFit;
                X(i,:) = GX(i,:);
            end
            if newFit < Silverback_Score
                Silverback_Score = newFit;
                Silverback = GX(i,:);
            end
        end

        %% Exploitation
        for i = 1:nPop
            if a >= w
                g = 2^C;
                delta = (abs(mean(GX)).^g).^(1/g);
                GX(i,:) = C*delta.*(X(i,:) - Silverback) + X(i,:);
            else
                if rand >= 0.5
                    h = randn(1, dim);
                else
                    h = randn(1,1);
                end
                r1 = rand;
                GX(i,:) = Silverback - (Silverback*(2*r1-1) - X(i,:)*(2*r1-1)) .* (Beta*h);
            end
        end

        GX = boundaryCheck(GX, lb, ub);

        % Group formation update
        for i = 1:nPop
            newFit = fobj(GX(i,:));
            if newFit < Pop_Fit(i)
                Pop_Fit(i) = newFit;
                X(i,:) = GX(i,:);
            end
            if newFit < Silverback_Score
                Silverback_Score = newFit;
                Silverback = GX(i,:);
            end
        end

        convergence_curve(It) = Silverback_Score;
    end

end

%% Helper function: Boundary check
function X = boundaryCheck(X, lb, ub)
    for i = 1:size(X,1)
        X(i,:) = max(min(X(i,:), ub), lb);
    end
end

%% Helper function: Population initialization
function X = initializePopulation(N, dim, lb, ub)
    if isscalar(lb)
        X = rand(N,dim).*(ub - lb) + lb;
    else
        X = zeros(N,dim);
        for i = 1:dim
            X(:,i) = rand(N,1).*(ub(i) - lb(i)) + lb(i);
        end
    end
end
