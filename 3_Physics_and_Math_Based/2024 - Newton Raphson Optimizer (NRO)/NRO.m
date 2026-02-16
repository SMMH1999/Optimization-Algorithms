function [bestFitness, bestPosition, convergenceCurve] = NRO(lb, ub, dim, nPop, maxItr, objFun)
    % ----------------------------------------------------------------------------------------------------------------------
    % Newton-Raphson-Based Optimizer (NRO)
    % Programmer: Dr. Pradeep Jangir
    % Reference: R. Sowmya, M. Premkumar, P. Jangir, Engineering Applications of Artificial Intelligence, Vol. 128, 107532, 2024
    % ----------------------------------------------------------------------------------------------------------------------
    % Description:
    % This function implements the Newton-Raphson-Based Optimizer (NRO) for continuous optimization problems.
    %
    % Input:
    %   lb       - Lower bound (scalar or 1 x dim vector)
    %   ub       - Upper bound (scalar or 1 x dim vector)
    %   dim      - Dimensionality of the search space
    %   nPop     - Population size
    %   maxItr   - Maximum number of iterations
    %   objFun   - Handle to objective function to minimize
    %
    % Output:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Best position vector found
    %   convergenceCurve - Convergence curve (best fitness per iteration)
    % ----------------------------------------------------------------------------------------------------------------------

    %% Parameters
    DF = 0.6; % Trap Avoidance deciding factor

    % Ensure bounds are vectors
    if isscalar(lb), lb = repmat(lb, 1, dim); end
    if isscalar(ub), ub = repmat(ub, 1, dim); end

    %% Initialization
    Position = initialization(nPop, dim, ub, lb);
    Fitness = zeros(nPop,1);
    for i = 1:nPop
        Fitness(i) = objFun(Position(i,:));
    end

    [FitnessSorted, Ind] = sort(Fitness);
    bestFitness = FitnessSorted(1);
    bestPosition = Position(Ind(1), :);
    Worst_Cost = FitnessSorted(end);
    Worst_Pos = Position(Ind(end), :);

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for it = 1:maxItr
        delta = (1 - (2*it/maxItr))^5; % Dynamic parameter

        for i = 1:nPop
            % Differential selection
            P = randperm(nPop,2); a1 = P(1); a2 = P(2);
            rho = rand * (bestPosition - Position(i,:)) + rand * (Position(a1,:) - Position(a2,:));

            % Newton-Raphson Search Rule
            Flag = 1;
            NRSR = SearchRule(bestPosition, Worst_Pos, Position(i,:), rho, Flag);
            X1 = Position(i,:) - NRSR + rho;
            X2 = bestPosition - NRSR + rho;

            % Update position
            Xupdate = zeros(1, dim);
            for j = 1:dim
                a1r = rand; a2r = rand;
                X3 = Position(i,j) - delta*(X2(j)-X1(j));
                Xupdate(j) = a1r*(a1r*X1(j) + (1-a2r)*X2(j)) + (1-a2r)*X3;
            end

            % Trap Avoidance Operator
            if rand < DF
                theta1 = -1 + 2*rand();
                theta2 = -0.5 + rand();
                beta = rand < 0.5;
                u1 = beta*3*rand + (1-beta);
                u2 = beta*rand + (1-beta);

                if u1 < 0.5
                    X_TAO = Xupdate + theta1*(u1*bestPosition - u2*Position(i,:)) + theta2*delta*(u1*mean(Position) - u2*Position(i,:));
                else
                    X_TAO = bestPosition + theta1*(u1*bestPosition - u2*Position(i,:)) + theta2*delta*(u1*mean(Position) - u2*Position(i,:));
                end
                Xnew = X_TAO;
            else
                Xnew = Xupdate;
            end

            % Boundary enforcement
            Xnew = min(max(Xnew, lb), ub);

            % Evaluate new solution
            Xnew_Cost = objFun(Xnew);

            % Update personal and global best
            if Xnew_Cost < Fitness(i)
                Position(i,:) = Xnew;
                Fitness(i) = Xnew_Cost;
                if Fitness(i) < bestFitness
                    bestPosition = Position(i,:);
                    bestFitness = Fitness(i);
                end
            end

            % Update global worst
            if Fitness(i) > Worst_Cost
                Worst_Pos = Position(i,:);
                Worst_Cost = Fitness(i);
            end
        end

        % Store convergence
        convergenceCurve(it) = bestFitness;
    end

end

%% Newton-Raphson Search Rule
function NRSR = SearchRule(Best_Pos, Worst_Pos, Position, rho, Flag)
    dim = numel(Position);
    DelX = rand(1, dim).*abs(Best_Pos - Position);

    % Initial Newton-Raphson step
    NRSR = randn*((Best_Pos - Worst_Pos).*DelX)./(2*(Best_Pos + Worst_Pos - 2*Position));

    if Flag == 1
        Xa = Position - NRSR + rho;
    else
        Xa = Best_Pos - NRSR + rho;
    end

    r1 = rand; r2 = rand;
    yp = r1*(mean(Xa + Position) + r1*DelX);
    yq = r2*(mean(Xa + Position) - r2*DelX);
    NRSR = randn*((yp - yq).*DelX)./(2*(yp + yq - 2*Position));
end

%% Population Initialization
function X = initialization(nP, dim, UB, LB)
    Boundary_no = numel(UB);
    if Boundary_no == 1
        X = rand(nP, dim).*(UB - LB) + LB;
    else
        X = zeros(nP, dim);
        for i = 1:dim
            X(:,i) = rand(nP,1).*(UB(i) - LB(i)) + LB(i);
        end
    end
end
