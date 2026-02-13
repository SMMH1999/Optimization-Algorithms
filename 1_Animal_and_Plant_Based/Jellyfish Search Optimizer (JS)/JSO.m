function [bestFitness, bestPosition, convergenceCurve] = JSO(lb, ub, dim, nPop, maxItr, objFun)
    %--------------------------------------------------------------------------%
    % Jellyfish Search Optimizer (JS)
    % Version: 1.0 | MATLAB R2016a+
    %
    % Authors:
    %   Professor Jui-Sheng Chou
    %   Ph.D. Candidate Dinh-Nhat Truong
    %
    % Reference:
    %   "A Novel Metaheuristic Optimizer Inspired By Behavior of Jellyfish
    %    in Ocean," Applied Mathematics and Computation, Volume 389, 15 Jan 2021, 125535.
    %    DOI: https://doi.org/10.1016/j.amc.2020.125535
    %
    % Description:
    %   Implementation of the Jellyfish Search Optimizer (JS) for global
    %   minimization problems. Supports scalar or vector bounds.
    %
    % Inputs:
    %   lb      - Lower bound (scalar or 1 x dim vector)
    %   ub      - Upper bound (scalar or 1 x dim vector)
    %   dim     - Number of decision variables
    %   nPop    - Population size
    %   maxItr  - Maximum number of iterations
    %   objFun  - Handle to objective function: fitness = objFun(position)
    %
    % Outputs:
    %   bestFitness       - Best fitness value found
    %   bestPosition      - Best solution vector
    %   convergenceCurve  - Best fitness at each iteration
    %--------------------------------------------------------------------------%

    %% Problem Definition
    VarSize = [1 dim];  % Decision variables size
    if isscalar(lb)
        VarMin = lb * ones(1, dim);
        VarMax = ub * ones(1, dim);
    else
        VarMin = lb;
        VarMax = ub;
    end

    %% Initialize Population (Logistic map)
    pop = initialization(nPop, dim, VarMax, VarMin);

    % Evaluate initial population
    popCost = zeros(nPop,1);
    for i = 1:nPop
        popCost(i) = objFun(pop(i,:), 0); % fnumber=0 placeholder
    end

    %% Main Loop
    convergenceCurve = zeros(1, maxItr);
    for it = 1:maxItr
        MeanPos = mean(pop,1);
        [sortedCost, idx] = sort(popCost);
        BestSol = pop(idx(1), :);
        BestCost = popCost(idx(1));

        for i = 1:nPop
            % Time control
            Ar = (1 - it / maxItr) * (2*rand - 1);
            if abs(Ar) >= 0.5
                % Following ocean current
                newSol = pop(i,:) + rand(VarSize).*(BestSol - 3*rand*MeanPos);
            else
                % Moving inside swarm
                if rand <= (1 - Ar)
                    j = i;
                    while j == i
                        j = randperm(nPop,1);
                    end
                    Step = pop(i,:) - pop(j,:);
                    if popCost(j) < popCost(i)
                        Step = -Step;
                    end
                    newSol = pop(i,:) + rand(VarSize).*Step;
                else
                    % Passive motion
                    newSol = pop(i,:) + 0.1*(VarMax-VarMin).*rand(VarSize);
                end
            end

            % Boundary check
            newSol = simplebounds(newSol, VarMin, VarMax);

            % Evaluation
            newCost = objFun(newSol, 0);
            if newCost < popCost(i)
                pop(i,:) = newSol;
                popCost(i) = newCost;
                if newCost < BestCost
                    BestCost = newCost;
                    BestSol = newSol;
                end
            end
        end

        % Record convergence
        convergenceCurve(it) = BestCost;

        % Early stopping
        if it >= 2000 && abs(convergenceCurve(it) - convergenceCurve(it-100)) < 1e-350
            convergenceCurve = convergenceCurve(1:it);
            break;
        end
    end

    bestPosition = BestSol;
    bestFitness = BestCost;

end

%% -------------------- Helper Functions -------------------- %%
function s = simplebounds(s, Lb, Ub)
    s = min(max(s, Lb), Ub);
end

function pop = initialization(num_pop, dim, Ub, Lb)
    % Logistic map initialization
    x = rand(1, dim);
    a = 4;
    pop = zeros(num_pop, dim);
    pop(1,:) = x;
    for i = 2:num_pop
        x = a*x.*(1-x);
        pop(i,:) = x;
    end
    % Scale to bounds
    for k = 1:dim
        pop(:,k) = Lb(k) + pop(:,k).*(Ub(k)-Lb(k));
    end
end
