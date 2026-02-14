function [bestFitness, bestPosition, convergenceCurve] = CDO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Chernobyl Disaster Optimizer (CDO)
    % =========================================================================
    % Author: H. Shehadeh (2023)
    % Refactored to benchmark single-function format
    %
    % Description:
    % Chernobyl Disaster Optimizer (CDO) is a population-based metaheuristic
    % inspired by radiation propagation after the Chernobyl disaster.
    % The algorithm models three leading radiations (Alpha, Beta, Gamma)
    % that guide the population toward promising regions in the search space.
    %
    % =========================================================================
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size (number of search agents)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % =========================================================================
    % Algorithm Parameters:
    %   a  - Linearly decreasing control parameter (from 3 to 0)
    %   a1 - Alpha radiation scaling (logarithmic random)
    %   a2 - Beta radiation scaling (logarithmic random)
    %   a3 - Gamma radiation scaling (logarithmic random)
    %
    % =========================================================================

    %% ========================= Initialization ===============================

    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    % Initialize Alpha, Beta, Gamma
    Alpha_pos = zeros(1, dim);
    Alpha_score = inf;

    Beta_pos = zeros(1, dim);
    Beta_score = inf;

    Gamma_pos = zeros(1, dim);
    Gamma_score = inf;

    % Initialize population
    Positions = rand(nPop, dim) .* (ub - lb) + lb;

    convergenceCurve = zeros(1, maxItr);

    %% =========================== Main Loop ==================================

    for t = 1:maxItr

        %% ----------- Fitness Evaluation & Leader Update ---------------------
        for i = 1:nPop

            % Boundary handling
            Positions(i,:) = max(Positions(i,:), lb);
            Positions(i,:) = min(Positions(i,:), ub);

            fitness = objFun(Positions(i,:));

            % Update Alpha
            if fitness < Alpha_score
                Alpha_score = fitness;
                Alpha_pos = Positions(i,:);
            end

            % Update Beta
            if fitness > Alpha_score && fitness < Beta_score
                Beta_score = fitness;
                Beta_pos = Positions(i,:);
            end

            % Update Gamma
            if fitness > Alpha_score && fitness > Beta_score && fitness < Gamma_score
                Gamma_score = fitness;
                Gamma_pos = Positions(i,:);
            end
        end

        %% ------------------ Control Parameters ------------------------------
        a = 3 - t * (3 / maxItr);

        a1 = log10((16000 - 1) * rand + 16000);
        a2 = log10((270000 - 1) * rand + 270000);
        a3 = log10((300000 - 1) * rand + 300000);

        %% ------------------ Position Update ---------------------------------
        for i = 1:nPop
            for j = 1:dim

                % ---------------- Alpha ----------------
                r1 = rand;
                r2 = rand;
                pa = pi * r1^2 / (0.25 * a1) - a * rand;
                C1 = r2^2 * pi;
                D_alpha = abs(C1 * Alpha_pos(j) - Positions(i,j));
                va = 0.25 * (Alpha_pos(j) - pa * D_alpha);

                % ---------------- Beta -----------------
                r1 = rand;
                r2 = rand;
                pb = pi * r1^2 / (0.5 * a2) - a * rand;
                C2 = r2^2 * pi;
                D_beta = abs(C2 * Beta_pos(j) - Positions(i,j));
                vb = 0.5 * (Beta_pos(j) - pb * D_beta);

                % ---------------- Gamma ----------------
                r1 = rand;
                r2 = rand;
                py = (pi * r1^2) / a3 - a * rand;
                C3 = r2^2 * pi;
                D_gamma = abs(C3 * Gamma_pos(j) - Positions(i,j));
                vy = Gamma_pos(j) - py * D_gamma;

                % Final update
                Positions(i,j) = (va + vb + vy) / 3;
            end
        end

        %% ------------------ Convergence Tracking -----------------------------
        convergenceCurve(t) = Alpha_score;

    end

    %% =========================== Output =====================================

    bestFitness  = Alpha_score;
    bestPosition = Alpha_pos;

end
