function [bestFitness, bestPosition, convergenceCurve] = Algorithm(lb, ub, dim, nPop, maxItr, objFun)
    %% Kirchhoff’s Law Algorithm (KLA)
    % Author: Mojtaba Ghasemi & Nima Khodadadi
    % Description: A physics-inspired, non-parametric metaheuristic algorithm
    %              for global optimization problems.
    %
    % INPUTS:
    %   lb       : scalar or vector, lower bounds of decision variables
    %   ub       : scalar or vector, upper bounds of decision variables
    %   dim      : integer, number of decision variables
    %   nPop     : integer, population size
    %   maxItr   : integer, maximum number of function evaluations
    %   objFun   : handle, objective function to minimize, objFun(x)
    %
    % OUTPUTS:
    %   bestFitness       : scalar, best cost value found
    %   bestPosition      : 1 x dim vector, position of best solution
    %   convergenceCurve  : maxItr x 1 vector, best fitness at each iteration

    %% Initialization
    if isscalar(lb)
        lb = repmat(lb, 1, dim);
    end
    if isscalar(ub)
        ub = repmat(ub, 1, dim);
    end

    ebs = realmin;  % small epsilon to avoid division by zero

    % Empty population template
    sol.Position = [];
    sol.Cost = [];

    % Initialize population
    pop = repmat(sol, nPop, 1);
    BestSol.Cost = inf;
    TT = -inf;

    for i = 1:nPop
        pop(i).Position = unifrnd(lb, ub, [1, dim]);
        pop(i).Cost = objFun(pop(i).Position);
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
        if pop(i).Cost > TT
            TT = pop(i).Cost;
        end
        BestCostHistory(i) = BestSol.Cost;
    end

    % Array to hold best cost per iteration
    convergenceCurve = zeros(maxItr, 1);

    %% Main Loop
    funcEvalCount = nPop;

    while funcEvalCount <= maxItr
        newpop = repmat(sol, nPop, 1);

        for i = 1:nPop
            newpop(i).Cost = inf;

            A = randperm(nPop);
            A(A == i) = [];
            a = A(1); b = A(2); jj = A(3);

            q  = ((pop(i).Cost - pop(jj).Cost) + ebs) / (abs(pop(i).Cost - pop(jj).Cost) + ebs);
            Q  = (pop(i).Cost - pop(a).Cost) / (abs(pop(i).Cost - pop(a).Cost) + ebs);
            Q2 = (pop(i).Cost - pop(b).Cost) / (abs(pop(i).Cost - pop(b).Cost) + ebs);

            q1  = (pop(jj).Cost / pop(i).Cost)^(2 * rand);
            Q1  = (pop(a).Cost / pop(i).Cost)^(2 * rand);
            Q21 = (pop(b).Cost / pop(i).Cost)^(2 * rand);

            S1 = q1 * q * rand([1, dim]) .* (pop(jj).Position - pop(i).Position);
            S2 = Q * Q1 * rand([1, dim]) .* (pop(a).Position - pop(i).Position);
            S3 = Q2 * Q21 * rand([1, dim]) .* (pop(b).Position - pop(i).Position);

            S = (rand + rand) * S1 + (rand + rand) * S2 + (rand + rand) * S3;

            newsol.Position = pop(i).Position + S;
            newsol.Position = max(newsol.Position, lb);
            newsol.Position = min(newsol.Position, ub);

            newsol.Cost = objFun(newsol.Position);

            if newsol.Cost <= pop(i).Cost
                pop(i) = newsol;
                if pop(i).Cost <= BestSol.Cost
                    BestSol = pop(i);
                end
            end

            funcEvalCount = funcEvalCount + 1;
            if funcEvalCount > maxItr
                break;
            end

            BestCostHistory(funcEvalCount) = BestSol.Cost;
        end

        convergenceCurve(funcEvalCount) = BestSol.Cost;
    end

    %% Output
    bestFitness = BestSol.Cost;
    bestPosition = BestSol.Position;

end
