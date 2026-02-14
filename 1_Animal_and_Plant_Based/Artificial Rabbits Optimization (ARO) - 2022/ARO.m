function [bestFitness, bestPosition, convergenceCurve] = ARO(lb, ub, dim, nPop, maxItr, objFun)
    %% Artificial Rabbits Optimization (ARO)
    % lb, ub : lower and upper bounds (scalar or vector)
    % dim    : problem dimensionality
    % nPop   : population size
    % maxItr : maximum number of iterations
    % objFun : handle to objective function
    %
    % Outputs:
    % bestFitness      : best fitness value found
    % bestPosition     : best position vector
    % convergenceCurve : history of best fitness values

    % Ensure lb and ub are vectors
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end

    % Initialize population
    PopPos = rand(nPop, dim) .* (ub - lb) + lb;
    PopFit = zeros(nPop, 1);
    for i = 1:nPop
        PopFit(i) = objFun(PopPos(i,:));
    end

    % Initialize global best
    [bestFitness, idx] = min(PopFit);
    bestPosition = PopPos(idx,:);

    convergenceCurve = zeros(maxItr,1);

    % Main loop
    for It = 1:maxItr
        theta = 2*(1 - It/maxItr);
        Direct1 = zeros(nPop, dim);
        Direct2 = zeros(nPop, dim);

        for i = 1:nPop
            % Running length L
            L = (exp(1) - exp(((It-1)/maxItr)^2)) * sin(2*pi*rand);

            % Direction vector 1
            rd = ceil(rand*dim);
            Direct1(i, randperm(dim, rd)) = 1;
            R = L .* Direct1(i,:);

            % Energy factor
            A = 2*log(1/rand) * theta;

            if A > 1
                % Exploration
                K = [1:i-1 i+1:nPop];
                RandInd = K(randi(nPop-1));
                newPos = PopPos(RandInd,:) + R .* (PopPos(i,:) - PopPos(RandInd,:)) ...
                    + round(0.5*(0.05+rand)) .* randn(1,dim);
            else
                % Exploitation / Hiding
                Direct2(i, ceil(rand*dim)) = 1;
                H = ((maxItr - It + 1)/maxItr) * randn;
                b = PopPos(i,:) + H * Direct2(i,:) .* PopPos(i,:);
                newPos = PopPos(i,:) + R .* (rand*b - PopPos(i,:));
            end

            % Apply bounds
            newPos = SpaceBound(newPos, ub, lb);

            % Evaluate
            newFit = objFun(newPos);

            % Update individual if improved
            if newFit < PopFit(i)
                PopFit(i) = newFit;
                PopPos(i,:) = newPos;
            end
        end

        % Update global best
        [minFit, idx] = min(PopFit);
        if minFit < bestFitness
            bestFitness = minFit;
            bestPosition = PopPos(idx,:);
        end

        % Store history
        convergenceCurve(It) = bestFitness;
    end
end

%% --- Helper: enforce bounds ---
function X = SpaceBound(X, Up, Low)
    Dim = length(X);
    S = (X > Up) + (X < Low);
    X = (rand(1,Dim).*(Up-Low)+Low).*S + X.*(~S);
end
