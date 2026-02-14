function [bestFitness, bestPosition, convergenceCurve] = ECO(lb, ub, dim, nPop, maxItr, objFun)
    %==========================================================================
    % Educational Competition Optimizer (ECO)
    % Authors: Junbo Lian, Ali Asghar Heidari, Huiling Chen
    %
    % Description:
    %   Population-based metaheuristic inspired by educational competition.
    %
    % Inputs:
    %   lb        - Lower bound of variables [scalar or 1 x dim]
    %   ub        - Upper bound of variables [scalar or 1 x dim]
    %   dim       - Dimension of the problem [scalar]
    %   nPop      - Population size [scalar]
    %   maxItr    - Maximum number of iterations [scalar]
    %   objFun    - Handle to objective function (fitness function)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Position of the best solution
    %   convergenceCurve - Best fitness at each iteration
    %==========================================================================

    %% Parameters
    H = 0.5;  % Learning habit boundary
    G1 = 0.2; % proportion of primary school
    G2 = 0.1; % proportion of middle school
    G1Number = round(nPop * G1);
    G2Number = round(nPop * G2);

    if isscalar(ub)
        ub = ub * ones(1, dim);
        lb = lb * ones(1, dim);
    end

    %% Initialization (Logistic chaotic map)
    X = initializationLogistic(nPop, dim, ub, lb);
    fitness = zeros(1, nPop);
    for i = 1:nPop
        fitness(i) = objFun(X(i,:));
    end
    [fitness, idx] = sort(fitness);
    X = X(idx, :);

    GBestF = fitness(1);
    GBestX = X(1, :);
    AveF = mean(fitness);

    convergenceCurve = zeros(1, maxItr);
    X_new = X;

    %% Main Loop
    for t = 1:maxItr
        R1 = rand(1);
        R2 = rand(1);
        P = 4 * randn * (1 - t / maxItr);
        E = (pi * t) / (P * maxItr);
        w = 0.1 * log(2 - (t / maxItr));

        for j = 1:nPop
            if mod(t,3) == 1
                if j <= G1Number
                    X_new(j,:) = X(j,:) + w * (mean(X(j,:)) - X(j,:)) .* Levy(dim);
                else
                    X_new(j,:) = X(j,:) + w * (close(X(j,:),1,X,G1Number) - X(j,:)) .* randn(1,dim);
                end
            elseif mod(t,3) == 2
                if j <= G2Number
                    X_new(j,:) = X(j,:) + (GBestX - mean(X)) * exp(t/maxItr - 1) .* Levy(dim);
                else
                    if R1 < H
                        X_new(j,:) = X(j,:) - w * close(X(j,:),2,X,G2Number) - P * (E*w*close(X(j,:),2,X,G2Number) - X(j,:));
                    else
                        X_new(j,:) = X(j,:) - w * close(X(j,:),2,X,G2Number) - P * (w*close(X(j,:),2,X,G2Number) - X(j,:));
                    end
                end
            else
                if j <= G2Number
                    X_new(j,:) = X(j,:) + (GBestX - X(j,:)) .* randn(1,dim) - (GBestX - X(j,:)) .* randn(1,dim);
                else
                    if R2 < H
                        X_new(j,:) = GBestX - P * (E * GBestX - X(j,:));
                    else
                        X_new(j,:) = GBestX - P * (GBestX - X(j,:));
                    end
                end
            end

            % Boundary control
            X_new(j,:) = max(min(X_new(j,:), ub), lb);

            % Evaluate fitness
            fitness_new(j) = objFun(X_new(j,:));
            if fitness_new(j) > fitness(j)
                fitness_new(j) = fitness(j);
                X_new(j,:) = X(j,:);
            end

            % Update global best
            if fitness_new(j) < GBestF
                GBestF = fitness_new(j);
                GBestX = X_new(j,:);
            end
        end

        X = X_new;
        fitness = fitness_new;
        convergenceCurve(t) = GBestF;

        % Sort for next iteration
        [fitness, idx] = sort(fitness);
        X = X(idx,:);
        AveF = mean(fitness);
    end

    bestFitness = GBestF;
    bestPosition = GBestX;
end
%% --- Sub-functions ---

% Levy flight
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    o = u ./ abs(v).^(1/beta);
end

% Choose nearest school
function o = close(t, G, X, NumG)
    m = X(1,:);
    for s = 1:NumG
        school = X(s,:);
        if sum(abs(m - t)) > sum(abs(school - t))
            m = school;
        end
    end
    o = m;
end

% Logistic chaotic map initialization
function Positions = initializationLogistic(pop, dim, ub, lb)
    Positions = zeros(pop, dim);
    for i = 1:pop
        for j = 1:dim
            x0 = rand;
            a = 4;
            x = a*x0*(1-x0);
            if isscalar(ub)
                Positions(i,j) = (ub - lb)*x + lb;
            else
                Positions(i,j) = (ub(j) - lb(j))*x + lb(j);
            end
            Positions(i,j) = min(max(Positions(i,j), lb(j)), ub(j));
            x0 = x;
        end
    end
end