function [bestFitness, bestPosition, convergenceCurve] = AOA(lb, ub, dim, nPop, maxItr, objFun)
    %% Arithmetic Optimization Algorithm (AOA)
    % Author: Original by Mirjalili et al.
    % Description:
    %   AOA is a metaheuristic optimization algorithm inspired by arithmetic operators.
    %
    % Inputs:
    %   lb      - Lower bound (scalar or 1xDim vector)
    %   ub      - Upper bound (scalar or 1xDim vector)
    %   dim     - Number of decision variables (scalar)
    %   nPop    - Population size (scalar)
    %   maxItr  - Maximum number of iterations (scalar)
    %   objFun  - Handle to the objective function, e.g., @(x) rastrigin(x)
    %
    % Outputs:
    %   bestFitness       - Best objective function value found
    %   bestPosition      - Position of the best solution (1xDim)
    %   convergenceCurve  - Best fitness at each iteration (1xMaxItr)

    %% Parameters
    MOP_Max = 1;
    MOP_Min = 0.2;
    Alpha = 5;
    Mu = 0.499;

    %% Initialization
    bestPosition = zeros(1, dim);
    bestFitness = inf;
    convergenceCurve = zeros(1, maxItr);

    X = initializePositions(nPop, dim, ub, lb);
    Xnew = X;

    F = zeros(1, nPop);
    Fnew = zeros(1, nPop);

    % Initial fitness evaluation
    for i = 1:nPop
        F(i) = objFun(X(i,:));
        if F(i) < bestFitness
            bestFitness = F(i);
            bestPosition = X(i,:);
        end
    end

    %% Main Loop
    for iter = 1:maxItr
        MOP = 1 - (iter^(1/Alpha) / maxItr^(1/Alpha));
        MOA = MOP_Min + iter * ((MOP_Max - MOP_Min)/maxItr);

        % Update positions
        for i = 1:nPop
            for j = 1:dim
                r1 = rand();
                if isscalar(lb) && isscalar(ub)
                    range = (ub - lb) * Mu + lb;
                else
                    range = (ub(j) - lb(j)) * Mu + lb(j);
                end

                if r1 < MOA
                    r2 = rand();
                    if r2 > 0.5
                        Xnew(i,j) = bestPosition(j) / (MOP+eps) * range;
                    else
                        Xnew(i,j) = bestPosition(j) * MOP * range;
                    end
                else
                    r3 = rand();
                    if r3 > 0.5
                        Xnew(i,j) = bestPosition(j) - MOP * range;
                    else
                        Xnew(i,j) = bestPosition(j) + MOP * range;
                    end
                end
            end

            % Boundary check
            if isscalar(lb) && isscalar(ub)
                Xnew(i,:) = max(min(Xnew(i,:), ub), lb);
            end

            % Fitness evaluation
            Fnew(i) = objFun(Xnew(i,:));
            if Fnew(i) < F(i)
                X(i,:) = Xnew(i,:);
                F(i) = Fnew(i);
            end
            if F(i) < bestFitness
                bestFitness = F(i);
                bestPosition = X(i,:);
            end
        end

        % Update convergence
        convergenceCurve(iter) = bestFitness;
    end

end

%% Local function: Initialization
function X = initializePositions(nPop, dim, ub, lb)
    if isscalar(ub) && isscalar(lb)
        X = rand(nPop, dim) * (ub - lb) + lb;
    else
        X = zeros(nPop, dim);
        for j = 1:dim
            X(:,j) = rand(nPop,1) * (ub(j) - lb(j)) + lb(j);
        end
    end
end
