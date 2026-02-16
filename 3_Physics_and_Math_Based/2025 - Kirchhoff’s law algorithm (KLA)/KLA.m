function [bestFitness, bestPosition, convergenceCurve] = KLA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Kirchhoff’s Law Algorithm (KLA)
    % =========================================================================
    % Author: Mojtaba Ghasemi & Nima Khodadadi
    %
    % Description:
    % Kirchhoff’s Law Algorithm (KLA) is a physics-inspired, non-parametric
    % metaheuristic optimization algorithm based on Kirchhoff’s electrical
    % circuit laws. Candidate solutions interact using voltage-like and
    % current-like relations derived from fitness differences to guide the
    % search toward promising regions of the search space.
    %
    % -------------------------------------------------------------------------
    % INPUTS:
    % lb        : (1×dim or scalar) Lower bound(s) of decision variables
    % ub        : (1×dim or scalar) Upper bound(s) of decision variables
    % dim       : (integer) Number of decision variables
    % nPop      : (integer) Population size
    % maxItr    : (integer) Maximum number of iterations
    % objFun    : (function handle) Objective function to minimize
    %
    % -------------------------------------------------------------------------
    % OUTPUTS:
    % bestFitness      : (scalar) Best objective value found
    % bestPosition     : (1×dim vector) Best solution found
    % convergenceCurve : (maxItr×1 vector) Best fitness at each iteration
    %
    % -------------------------------------------------------------------------
    % Notes:
    % - Minimization problem is assumed.
    % - Supports scalar or vector bounds.
    % - Boundary constraints are strictly enforced.
    % - No internal display or printing.
    % =========================================================================

    %% ========================== Initialization ==============================
    if isscalar(lb)
        lb = repmat(lb, 1, dim);
    end
    if isscalar(ub)
        ub = repmat(ub, 1, dim);
    end

    ebs = realmin;  % Small epsilon to avoid division by zero

    % Population structure
    emptySol.Position = [];
    emptySol.Cost     = [];

    pop = repmat(emptySol, nPop, 1);
    BestSol.Cost = inf;

    % Initialize population
    for i = 1:nPop
        pop(i).Position = lb + rand(1, dim) .* (ub - lb);
        pop(i).Cost     = objFun(pop(i).Position);

        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end

    convergenceCurve = zeros(maxItr, 1);

    %% ============================ Main Loop =================================
    for t = 1:maxItr

        for i = 1:nPop

            % Randomly select three different individuals
            idx = randperm(nPop);
            idx(idx == i) = [];
            a  = idx(1);
            b  = idx(2);
            jj = idx(3);

            % Kirchhoff-inspired coefficients
            q  = ((pop(i).Cost - pop(jj).Cost) + ebs) / ...
                (abs(pop(i).Cost - pop(jj).Cost) + ebs);

            Q  = (pop(i).Cost - pop(a).Cost) / ...
                (abs(pop(i).Cost - pop(a).Cost) + ebs);

            Q2 = (pop(i).Cost - pop(b).Cost) / ...
                (abs(pop(i).Cost - pop(b).Cost) + ebs);

            q1  = (pop(jj).Cost / pop(i).Cost)^(2 * rand);
            Q1  = (pop(a).Cost  / pop(i).Cost)^(2 * rand);
            Q21 = (pop(b).Cost  / pop(i).Cost)^(2 * rand);

            % Step components
            S1 = q1  * q  * rand(1, dim) .* (pop(jj).Position - pop(i).Position);
            S2 = Q   * Q1 * rand(1, dim) .* (pop(a).Position  - pop(i).Position);
            S3 = Q2  * Q21* rand(1, dim) .* (pop(b).Position  - pop(i).Position);

            S = (rand + rand) * S1 + ...
                (rand + rand) * S2 + ...
                (rand + rand) * S3;

            % New candidate solution
            newPos = pop(i).Position + S;

            % Boundary control
            newPos = max(newPos, lb);
            newPos = min(newPos, ub);

            newCost = objFun(newPos);

            % Greedy selection
            if newCost <= pop(i).Cost
                pop(i).Position = newPos;
                pop(i).Cost     = newCost;

                if newCost <= BestSol.Cost
                    BestSol = pop(i);
                end
            end
        end

        % Store best fitness of this iteration
        convergenceCurve(t) = BestSol.Cost;
    end

    %% ============================== Output ==================================
    bestFitness  = BestSol.Cost;
    bestPosition = BestSol.Position;

end
