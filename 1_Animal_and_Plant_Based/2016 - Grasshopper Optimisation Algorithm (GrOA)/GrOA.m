function [bestFitness, bestPosition, convergenceCurve] = GOA(lb, ub, dim, nPop, maxItr, objFun)
    %_________________________________________________________________________%
    %  Grasshopper Optimization Algorithm (GOA) - Unified Benchmark Version    %
    %                                                                         %
    %  Author & Programmer: Seyedali Mirjalili                                 %
    %  Original Paper: S. Saremi, S. Mirjalili, A. Lewis,                     %
    %                  "Grasshopper Optimisation Algorithm: Theory and        %
    %                   Application", Advances in Engineering Software, 2017  %
    %  DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004               %
    %_________________________________________________________________________%
    %
    % DESCRIPTION:
    %   This function implements the Grasshopper Optimization Algorithm (GOA)
    %   for continuous optimization problems.
    %
    % INPUTS:
    %   lb      - Lower bound (scalar or vector) [1 x dim] or scalar
    %   ub      - Upper bound (scalar or vector) [1 x dim] or scalar
    %   dim     - Number of dimensions (problem variables)
    %   nPop    - Population size (number of grasshoppers)
    %   maxItr  - Maximum number of iterations
    %   objFun  - Function handle of the objective function to minimize
    %
    % OUTPUTS:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Position of the best solution
    %   convergenceCurve - Array of bestFitness values at each iteration
    %
    % TUNABLE PARAMETERS:
    %   cMax = 1, cMin = 0.00004 : control the intensity of social interaction

    %% Initialization
    flag = 0;  % flag to handle odd dimensions
    if isscalar(ub)
        ub = ones(dim,1)*ub;
        lb = ones(dim,1)*lb;
    end

    if rem(dim,2) ~= 0
        dim = dim + 1;
        ub = [ub; 100];
        lb = [lb; -100];
        flag = 1;
    end

    % Initialize population
    Positions = initialization(nPop, dim, ub, lb);
    Fitness = zeros(1, nPop);

    convergenceCurve = zeros(1, maxItr);

    % Control parameters
    cMax = 1;
    cMin = 0.00004;

    %% Evaluate initial fitness
    for i = 1:nPop
        if flag
            Fitness(i) = objFun(Positions(i,1:end-1));
        else
            Fitness(i) = objFun(Positions(i,:));
        end
    end

    [sortedFitness, idx] = sort(Fitness);
    SortedPositions = Positions(idx, :);

    bestPosition = SortedPositions(1,:);
    bestFitness = sortedFitness(1);

    %% Main loop
    for iter = 2:maxItr
        c = cMax - iter*((cMax - cMin)/maxItr); % social force coefficient

        for i = 1:nPop
            S_total = zeros(dim,1);
            for k = 1:2:dim
                S_i = zeros(2,1);
                for j = 1:nPop
                    if i ~= j
                        Dist = distance(Positions(i,k:k+1), Positions(j,k:k+1));
                        r_ij_vec = (Positions(j,k:k+1) - Positions(i,k:k+1))/(Dist+eps);
                        xj_xi = 2 + rem(Dist,2);
                        s_ij = ((ub(k:k+1)-lb(k:k+1))*c/2) .* S_func(xj_xi) .* r_ij_vec;
                        S_i = S_i + s_ij;
                    end
                end
                S_total(k:k+1) = S_i;
            end
            X_new = c*S_total' + bestPosition;
            Positions(i,:) = X_new;
        end

        % Boundary control and fitness evaluation
        for i = 1:nPop
            Positions(i,:) = min(max(Positions(i,:), lb'), ub');

            if flag
                Fitness(i) = objFun(Positions(i,1:end-1));
            else
                Fitness(i) = objFun(Positions(i,:));
            end

            % Update best solution
            if Fitness(i) < bestFitness
                bestFitness = Fitness(i);
                bestPosition = Positions(i,:);
            end
        end

        convergenceCurve(iter) = bestFitness;
    end

    if flag
        bestPosition = bestPosition(1:end-1);
    end

end

%% ----------------- Helper Functions ----------------- %%
function X = initialization(N, dim, ub, lb)
    if size(ub,1) == 1
        X = rand(N, dim) .* (ub - lb) + lb;
    else
        X = zeros(N, dim);
        for i = 1:dim
            high = ub(i); low = lb(i);
            X(:,i) = rand(N,1) .* (high-low) + low;
        end
    end
end

function d = distance(a,b)
    d = sqrt(sum((a-b).^2));
end

function o = S_func(r)
    f = 0.5;
    l = 1.5;
    o = f*exp(-r/l) - exp(-r);
end
