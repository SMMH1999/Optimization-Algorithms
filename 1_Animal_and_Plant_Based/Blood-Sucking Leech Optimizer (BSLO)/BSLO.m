function [bestFitness, bestPosition, convergenceCurve] = BSLO(lb, ub, dim, nPop, maxItr, objFun)
    %% Blood-Sucking Leech Optimizer (BSLO)
    % Abbreviation: BSLO
    %
    % Author: Jianfu Bai
    % Email: Jianfu.Bai@UGent.be, magd.abdelwahab@ugent.be
    % Soete Laboratory, Ghent University, Belgium
    %
    % Reference: Jianfu Bai, H. Nguyen-Xuan, Elena Atroshchenko, Gregor Kosec,
    % Lihua Wang, Magd Abdel Wahab, "Blood-sucking leech optimizer",
    % Advances in Engineering Software, 2024, 195:103696.
    % DOI: 10.1016/j.advengsoft.2024.103696
    %
    % Description:
    %   Implements the Blood-Sucking Leech Optimizer (BSLO) for continuous
    %   minimization problems.
    %
    % Inputs:
    %   lb      - Lower bound (1xD vector or scalar)
    %   ub      - Upper bound (1xD vector or scalar)
    %   dim     - Number of decision variables (integer)
    %   nPop    - Number of search agents (integer)
    %   maxItr  - Maximum number of iterations (integer)
    %   objFun  - Objective function handle (@f)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found (scalar)
    %   bestPosition     - Position of the best solution (1xD vector)
    %   convergenceCurve - Best fitness at each iteration (1xmaxItr vector)
    %
    % Tunable parameters:
    %   m  - Fraction of directional leeches (default: 0.8)
    %   a  - Probability factor for exploitation/exploration (default: 0.97)
    %   b  - Small weight for exploration (default: 0.001)
    %   t1 - Re-tracking start iteration (default: 20)
    %   t2 - Re-tracking look-back period (default: 20)

    %% -------------------- Initialization --------------------
    bestPosition = zeros(1, dim);
    bestFitness = inf;
    convergenceCurve = zeros(1, maxItr);
    Temp_best_fitness = zeros(1, maxItr);

    % Initialize population
    Positions = initializePopulation(nPop, dim, ub, lb);
    fitness = zeros(1, nPop);

    % Algorithm parameters
    t = 0; m = 0.8; a = 0.97; b = 0.001; t1 = 20; t2 = 20;

    %% -------------------- Main Loop --------------------
    while t < maxItr
        N1 = floor((m + (1 - m) * (t / maxItr)^2) * nPop);

        % Evaluate fitness and update best
        for i = 1:nPop
            % Boundary check
            Positions(i,:) = max(min(Positions(i,:), ub), lb);

            % Fitness evaluation
            fitness(i) = objFun(Positions(i,:));

            % Update best solution
            if fitness(i) <= bestFitness
                bestFitness = fitness(i);
                bestPosition = Positions(i,:);
            end
        end

        Prey_Position = bestPosition;

        % Re-tracking strategy
        Temp_best_fitness(t+1) = bestFitness;
        if t > t1 && Temp_best_fitness(t+1) == Temp_best_fitness(t+1-t2)
            for i = 1:nPop
                if fitness(i) == bestFitness
                    Positions(i,:) = rand(1, dim) .* (ub - lb) + lb;
                end
            end
        end

        % Dynamic coefficients
        if rand() < 0.5
            s = 8 - 1 * (-(t / maxItr)^2 + 1);
        else
            s = 8 - 7 * (-(t / maxItr)^2 + 1);
        end
        beta = -0.5 * (t / maxItr)^6 + (t / maxItr)^4 + 1.5;
        LV = 0.5 * levy(nPop, dim, beta);

        % Random integer matrices
        minValue = 1;
        maxValue = floor(nPop * (1 + t / maxItr));
        k2 = randi([minValue, maxValue], nPop, dim);
        k = randi([minValue, dim], nPop, dim);

        %% Directional Leeches Update
        for i = 1:N1
            for j = 1:dim
                r1 = 2 * rand() - 1;
                r2 = 2 * rand() - 1;
                r3 = 2 * rand() - 1;

                PD = s * (1 - (t / maxItr)) * r1;

                if abs(PD) >= 1
                    % Exploration
                    W1 = (1 - t / maxItr) * b * LV(i,j);
                    L1 = r2 * abs(Prey_Position(j) - Positions(i,j)) * PD * (1 - k2(i,j)/nPop);
                    L2 = abs(Prey_Position(j) - Positions(i,k(i,j))) * PD * (1 - (r2^2) * (k2(i,j)/nPop));

                    if rand() < a
                        if abs(Prey_Position(j)) > abs(Positions(i,j))
                            Positions(i,j) = Positions(i,j) + W1 * Positions(i,j) - L1;
                        else
                            Positions(i,j) = Positions(i,j) + W1 * Positions(i,j) + L1;
                        end
                    else
                        if abs(Prey_Position(j)) > abs(Positions(i,j))
                            Positions(i,j) = Positions(i,j) + W1 * Positions(i,k(i,j)) - L2;
                        else
                            Positions(i,j) = Positions(i,j) + W1 * Positions(i,k(i,j)) + L2;
                        end
                    end
                else
                    % Exploitation
                    if t >= 0.1 * maxItr
                        b = 0.00001;
                    end
                    W1 = (1 - t / maxItr) * b * LV(i,j);
                    L3 = abs(Prey_Position(j) - Positions(i,j)) * PD * (1 - r3 * k2(i,j)/nPop);
                    L4 = abs(Prey_Position(j) - Positions(i,k(i,j))) * PD * (1 - r3 * k2(i,j)/nPop);

                    if rand() < a
                        if abs(Prey_Position(j)) > abs(Positions(i,j))
                            Positions(i,j) = Prey_Position(j) + W1 * Prey_Position(j) - L3;
                        else
                            Positions(i,j) = Prey_Position(j) + W1 * Prey_Position(j) + L3;
                        end
                    else
                        if abs(Prey_Position(j)) > abs(Positions(i,j))
                            Positions(i,j) = Prey_Position(j) + W1 * Prey_Position(j) - L4;
                        else
                            Positions(i,j) = Prey_Position(j) + W1 * Prey_Position(j) + L4;
                        end
                    end
                end
            end
        end

        %% Directionless Leeches Update
        for i = N1+1:nPop
            for j = 1:dim
                if min(lb) >= 0
                    LV(i,j) = abs(LV(i,j));
                end
                if rand() > 0.5
                    Positions(i,j) = (t / maxItr) * LV(i,j) * Positions(i,j) * abs(Prey_Position(j) - Positions(i,j));
                else
                    Positions(i,j) = (t / maxItr) * LV(i,j) * Prey_Position(j) * abs(Prey_Position(j) - Positions(i,j));
                end
            end
        end

        t = t + 1;
        convergenceCurve(t) = bestFitness;
    end

    %% -------------------- Helper Functions --------------------
    function X = initializePopulation(nAgents, dim, ub, lb)
        if numel(ub) == 1
            X = rand(nAgents, dim) * (ub - lb) + lb;
        else
            X = zeros(nAgents, dim);
            for d = 1:dim
                X(:,d) = rand(nAgents,1) * (ub(d) - lb(d)) + lb(d);
            end
        end
    end

    function Z = levy(n, m, beta)
        % Generate Levy flight steps
        num = gamma(1 + beta) * sin(pi * beta / 2);
        den = gamma((1 + beta)/2) * beta * 2^((beta - 1)/2);
        sigma_u = (num / den)^(1/beta);
        u = random('Normal', 0, sigma_u, n, m);
        v = random('Normal', 0, 1, n, m);
        Z = u ./ (abs(v).^(1/beta));
    end

end
