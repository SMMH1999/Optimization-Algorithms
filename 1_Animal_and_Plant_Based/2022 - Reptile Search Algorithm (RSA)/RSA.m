function [bestFitness, bestPosition, convergenceCurve] = RSA(lb, ub, dim, nPop, maxItr, objFun)
    %% Reptile Search Algorithm (RSA)
    % Abbreviation: RSA
    %
    % Author: Laith Abualigah
    % Main reference: "Reptile Search Algorithm (RSA): A novel nature-inspired metaheuristic algorithm"
    %
    % Description:
    %   This function performs global optimization using the Reptile Search Algorithm (RSA),
    %   a nature-inspired metaheuristic mimicking the hunting behavior of reptiles.
    %
    % Inputs:
    %   lb      - scalar or 1 x dim vector, lower bounds of decision variables
    %   ub      - scalar or 1 x dim vector, upper bounds of decision variables
    %   dim     - integer, dimension of the problem
    %   nPop    - integer, population size
    %   maxItr  - integer, maximum number of iterations
    %   objFun  - function handle, objective function to minimize
    %
    % Outputs:
    %   bestFitness       - best objective function value found
    %   bestPosition      - 1 x dim vector of decision variables corresponding to bestFitness
    %   convergenceCurve  - 1 x maxItr vector storing the bestFitness at each iteration
    %
    % Tunable parameters:
    %   Alpha  - controls exploitation step (default 0.1)
    %   Beta   - controls exploration step (default 0.005)

    %% Initialization
    bestPosition = zeros(1, dim);
    bestFitness = inf;

    % Initialize population
    X = initializePopulation(nPop, dim, lb, ub);
    Xnew = zeros(nPop, dim);

    % Convergence curve
    convergenceCurve = zeros(1, maxItr);

    % Fitness arrays
    F = zeros(1, nPop);

    % Parameters
    Alpha = 0.1;
    Beta = 0.005;

    % Evaluate initial population
    for i = 1:nPop
        F(i) = objFun(X(i,:));
        if F(i) < bestFitness
            bestFitness = F(i);
            bestPosition = X(i,:);
        end
    end

    %% Main loop
    for t = 1:maxItr
        ES = 2*randn*(1 - t/maxItr); % Probability ratio
        for i = 2:nPop
            for j = 1:dim
                R   = bestPosition(j) - X(randi([1 nPop]), j) / (bestPosition(j) + eps);
                P   = Alpha + (X(i,j) - mean(X(i,:))) / (bestPosition(j)*(ub(j)-lb(j)) + eps);
                Eta = bestPosition(j) * P;

                % Position update based on iteration stage
                if t < maxItr/4
                    Xnew(i,j) = bestPosition(j) - Eta*Beta - R*rand;
                elseif t < maxItr/2
                    Xnew(i,j) = bestPosition(j) * X(randi([1 nPop]), j) * ES * rand;
                elseif t < 3*maxItr/4
                    Xnew(i,j) = bestPosition(j) * P * rand;
                else
                    Xnew(i,j) = bestPosition(j) - Eta*eps - R*rand;
                end
            end

            % Boundary check
            Xnew(i,:) = max(min(Xnew(i,:), ub), lb);

            % Evaluate new solution
            Fnew = objFun(Xnew(i,:));
            if Fnew < F(i)
                X(i,:) = Xnew(i,:);
                F(i) = Fnew;
            end

            % Update global best
            if F(i) < bestFitness
                bestFitness = F(i);
                bestPosition = X(i,:);
            end
        end

        % Store convergence
        convergenceCurve(t) = bestFitness;
    end
end

%% Population initialization function
function X = initializePopulation(N, Dim, UB, LB)
    % Supports scalar or vector bounds
    if isscalar(UB)
        X = rand(N, Dim) .* (UB - LB) + LB;
    else
        X = zeros(N, Dim);
        for i = 1:Dim
            X(:,i) = rand(N,1) .* (UB(i)-LB(i)) + LB(i);
        end
    end
end
