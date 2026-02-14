function [bestFitness, bestPosition, convergenceCurve] = ES(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Evolution Strategy (ES)
    % Classical (mu,lambda) / (mu+lambda) Evolution Strategy based on:
    % Thomas Back, "Evolutionary Algorithms in Theory and Practice", 1996.
    %
    % Refactored and unified into benchmark format.
    %
    % INPUTS:
    %   lb        : Lower bound (scalar or 1×dim vector)
    %   ub        : Upper bound (scalar or 1×dim vector)
    %   dim       : Problem dimension (number of decision variables)
    %   nPop      : Parent population size (mu)
    %   maxItr    : Maximum number of generations
    %   objFun    : Objective function handle (minimization)
    %
    % OUTPUTS:
    %   bestFitness      : Best fitness value found (scalar)
    %   bestPosition     : Best decision vector found (1×dim)
    %   convergenceCurve : Best fitness at each iteration (maxItr×1)
    %
    % INTERNAL PARAMETERS (fixed to preserve classical ES behavior):
    %   lambda  = 7*mu          : Offspring size
    %   sel     = '+'           : Selection scheme ('+' or ',')
    %   rec_obj = 4             : Intermediate recombination (object variables)
    %   rec_str = 4             : Intermediate recombination (strategy params)
    %
    % NOTES:
    %   - Fully self-adaptive covariance matrix ES.
    %   - Supports scalar/vector bounds.
    %   - Minimization problem.
    % =========================================================================

    %% -------------------- Parameter Definition -----------------------------
    mu      = nPop;
    lambda  = 7 * mu;      % standard ES choice
    sel     = '+';         % (mu + lambda) selection
    rec_obj = 4;           % intermediate recombination
    rec_str = 4;

    %% -------------------- Bound Handling -----------------------------------
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end
    lb = lb(:)';
    ub = ub(:)';

    %% -------------------- Initialization ------------------------------------
    X = zeros(dim, mu);
    for i = 1:dim
        X(i,:) = lb(i) + rand(1,mu).*(ub(i)-lb(i));
    end

    sigma = cell(1,mu);
    for i = 1:mu
        tmp = rand(dim);
        sigma{i} = tmp + tmp';
    end

    alpha = cell(1,mu);
    for i = 1:mu
        alpha{i} = zeros(dim);
        for j = 1:dim-1
            for k = j+1:dim
                alpha{i}(j,k) = 0.5 * atan2(2*sigma{i}(j,k), ...
                    (sigma{i}(j,j)^2 - sigma{i}(k,k)^2));
            end
        end
    end

    fitness = zeros(1,mu);
    for i = 1:mu
        fitness(i) = objFun(X(:,i)');
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = X(:,idx)';
    convergenceCurve = zeros(maxItr,1);
    convergenceCurve(1) = bestFitness;

    %% -------------------- Main Loop ----------------------------------------
    for t = 2:maxItr

        %% -------- Recombination --------------------------------------------
        Xr = zeros(dim,lambda);
        sigmar = cell(1,lambda);
        alphar = cell(1,lambda);

        for i = 1:lambda
            parents = randperm(mu,2);
            % intermediate recombination (object)
            Xr(:,i) = X(:,parents(1)) + ...
                0.5*(X(:,parents(2)) - X(:,parents(1)));

            % intermediate recombination (strategy)
            sigmar{i} = sigma{parents(1)} + ...
                0.5*(sigma{parents(2)} - sigma{parents(1)});
            alphar{i} = alpha{parents(1)} + ...
                0.5*(alpha{parents(2)} - alpha{parents(1)});
        end

        %% -------- Mutation --------------------------------------------------
        tau  = 1/(sqrt(2*sqrt(dim)));
        taup = 1/(sqrt(2*dim));
        beta = 5*pi/180;

        Xm = zeros(dim,lambda);
        sigmam = cell(1,lambda);
        alpham = cell(1,lambda);

        for i = 1:lambda
            tmp = randn(dim);
            sigmam{i} = sigmar{i} .* exp(taup*randn + tau*(tmp+tmp'));

            tmp2 = rand(dim);
            alpham{i} = alphar{i} + beta*triu((tmp2+tmp2'),1);

            R = eye(dim);
            for m = 1:dim-1
                for q = m+1:dim
                    T = eye(dim);
                    T([m q],[m q]) = ...
                        [cos(alpham{i}(m,m)) -sin(alpham{i}(m,q));
                        sin(alpham{i}(q,m))  cos(alpham{i}(q,q))];
                    R = R*T;
                end
            end

            Xm(:,i) = Xr(:,i) + ...
                R*sqrt(diag(diag(sigmam{i}))) * randn(dim,1);

            Xm(:,i) = max(Xm(:,i), lb');
            Xm(:,i) = min(Xm(:,i), ub');
        end

        %% -------- Evaluation -----------------------------------------------
        fitOff = zeros(1,lambda);
        for i = 1:lambda
            fitOff(i) = objFun(Xm(:,i)');
        end

        %% -------- Selection -------------------------------------------------
        if sel == ','      % (mu,lambda)
            [~, idx] = sort(fitOff);
            X = Xm(:,idx(1:mu));
            sigma = sigmam(idx(1:mu));
            alpha = alpham(idx(1:mu));
            fitness = fitOff(idx(1:mu));
        else               % (mu+lambda)
            Xaug = [Xm X];
            sigmaAug = [sigmam sigma];
            alphaAug = [alpham alpha];
            fitAug = [fitOff fitness];
            [~, idx] = sort(fitAug);
            X = Xaug(:,idx(1:mu));
            sigma = sigmaAug(idx(1:mu));
            alpha = alphaAug(idx(1:mu));
            fitness = fitAug(idx(1:mu));
        end

        %% -------- Best Update ----------------------------------------------
        [currentBest, idx] = min(fitness);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(:,idx)';
        end

        convergenceCurve(t) = bestFitness;
    end

end
