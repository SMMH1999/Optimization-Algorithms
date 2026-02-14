function [bestFitness, bestPosition, convergenceCurve] = PLO(lb, ub, dim, nPop, maxItr, objFun)
    %% Polar Lights Optimizer (PLO)
    % Author: Chong Yuan, Dong Zhao, Ali Asghar Heidari, Lei Liu, Yi Chen, Huiling Chen
    % Description: PLO is a bio-inspired optimizer mimicking polar light behavior.
    %
    % Inputs:
    %   lb        : lower bound (scalar or 1xdim vector)
    %   ub        : upper bound (scalar or 1xdim vector)
    %   dim       : number of decision variables
    %   nPop      : number of search agents (population size)
    %   maxItr    : maximum number of iterations
    %   objFun    : handle to objective function
    %
    % Outputs:
    %   bestFitness       : best objective value found
    %   bestPosition      : position of the best solution
    %   convergenceCurve  : best fitness value at each iteration

    %% Initialization
    X = initialization(nPop, dim, ub, lb);   % population positions
    V = ones(nPop, dim);                     % velocity factor
    fitness = inf(nPop, 1);

    for i = 1:nPop
        fitness(i) = objFun(X(i,:));
    end

    [fitness, SortOrder] = sort(fitness);
    X = X(SortOrder, :);
    bestPosition = X(1,:);
    bestFitness = fitness(1);

    convergenceCurve = bestFitness;

    %% Main loop
    for it = 2:maxItr

        X_sum = sum(X, 1);
        X_mean = X_sum / nPop;
        w1 = tansig((it/maxItr)^4);
        w2 = exp(-(2*it/maxItr)^3);

        X_new = zeros(nPop, dim);

        for i = 1:nPop
            a = rand()/2 + 1;
            V(i,:) = exp((1 - a)/100 * it);
            LS = V(i,:);
            GS = Levy(dim) .* (X_mean - X(i,:) + (lb + rand(1,dim).*(ub - lb))/2);
            X_new(i,:) = X(i,:) + (w1*LS + w2*GS) .* rand(1,dim);
        end

        E = sqrt(it / maxItr);
        A = randperm(nPop);

        for i = 1:nPop
            for j = 1:dim
                if rand < 0.05 && rand < E
                    X_new(i,j) = X(i,j) + sin(rand*pi)*(X(i,j) - X(A(i), j));
                end
            end
            % Boundary control
            X_new(i,:) = max(min(X_new(i,:), ub), lb);

            f_new = objFun(X_new(i,:));
            if f_new < fitness(i)
                X(i,:) = X_new(i,:);
                fitness(i) = f_new;
            end
        end

        [fitness, SortOrder] = sort(fitness);
        X = X(SortOrder, :);
        if fitness(1) < bestFitness
            bestPosition = X(1,:);
            bestFitness = fitness(1);
        end

        convergenceCurve(it) = bestFitness;
    end

end

%% Levy flight generator
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    o = u ./ abs(v).^(1/beta);
end

%% Population initialization
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Boundary_no = numel(ub);

    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    else
        Positions = zeros(SearchAgents_no, dim);
        for i = 1:dim
            Positions(:,i) = rand(SearchAgents_no,1) .* (ub(i) - lb(i)) + lb(i);
        end
    end
end
