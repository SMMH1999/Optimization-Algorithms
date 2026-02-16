function [bestFitness, bestPosition, convergenceCurve] = STOA(lb, ub, dim, nPop, maxItr, objFun)
    %################################################################################################%
    % STOA: Seagull-Tracking Optimization Algorithm (STOA)                                           %
    % Authors: Gaurav Dhiman and Amandeep Kaur                                                      %
    % DOI: https://doi.org/10.1016/j.engappai.2019.03.021                                           %
    % Link: https://www.sciencedirect.com/science/article/abs/pii/S0952197619300715                 %
    % Description:                                                                                  %
    %   STOA is a bio-inspired optimization algorithm designed for industrial engineering problems. %
    %   It simulates the hunting strategy and movement patterns of seagulls.                        %
    %                                                                                                %
    % Inputs:                                                                                        %
    %   lb        - Lower bound (scalar or 1 x dim vector)                                         %
    %   ub        - Upper bound (scalar or 1 x dim vector)                                         %
    %   dim       - Dimension of the problem (integer)                                             %
    %   nPop      - Number of search agents (population size)                                      %
    %   maxItr    - Maximum number of iterations                                                   %
    %   objFun    - Handle to objective function, objFun(x)                                        %
    %                                                                                                %
    % Outputs:                                                                                       %
    %   bestFitness      - Best objective function value found                                      %
    %   bestPosition     - Position vector corresponding to bestFitness                             %
    %   convergenceCurve - Array containing bestFitness at each iteration                            %
    %                                                                                                %
    % Tunable Parameters:                                                                            %
    %   Sa - Seagull attack coefficient, linearly decreasing from 2 to 0                              %
    %   b  - Spiral shape coefficient for movement update, fixed at 1                               %
    %################################################################################################%

    %---------------------------- Initialization ----------------------------%
    bestPosition = zeros(1, dim);
    bestFitness = inf;

    positions = initializePopulation(nPop, dim, ub, lb);

    convergenceCurve = zeros(1, maxItr);

    %---------------------------- Main Optimization Loop --------------------%
    for iter = 1:maxItr

        %-------------------- Boundary Handling & Fitness Evaluation --------------------%
        for i = 1:nPop
            positions(i,:) = max(min(positions(i,:), ub), lb);  % enforce bounds
            fitness = objFun(positions(i,:));
            if fitness < bestFitness
                bestFitness = fitness;
                bestPosition = positions(i,:);
            end
        end

        %-------------------- Seagull Attack Coefficient --------------------%
        Sa = 2 - iter * (2 / maxItr);

        %-------------------- Position Update (Spiral + Attack) --------------------%
        for i = 1:nPop
            for j = 1:dim
                r1 = 0.5 * rand();
                r2 = 0.5 * rand();
                X1 = 2 * Sa * r1 - Sa;
                b = 1;
                ll = (Sa - 1) * rand() + 1;
                D_alpha = Sa * positions(i,j) + X1 * (bestPosition(j) - positions(i,j));
                positions(i,j) = D_alpha * exp(b * ll) * cos(2 * pi * ll) + bestPosition(j);
            end
        end

        %-------------------- Store Convergence --------------------%
        convergenceCurve(iter) = bestFitness;
    end

end

%==================== Helper Function: Initialize Population ====================%
function pos = initializePopulation(nPop, dim, ub, lb)
    if isscalar(ub) && isscalar(lb)
        pos = rand(nPop, dim) * (ub - lb) + lb;
    else
        pos = zeros(nPop, dim);
        for d = 1:dim
            pos(:,d) = rand(nPop,1) * (ub(d) - lb(d)) + lb(d);
        end
    end
end
