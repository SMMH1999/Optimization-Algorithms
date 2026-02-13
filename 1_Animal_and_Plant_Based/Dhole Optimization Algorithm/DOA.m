function [bestFitness, bestPosition, convergenceCurve] = DOA(lb, ub, dim, nPop, maxItr, objFun)
    % -------------------------------------------------------------------------
    % Dhole Optimization Algorithm (DOA)
    % -------------------------------------------------------------------------
    % Author: (Original version not specified in provided code)
    % Refactored to Benchmark / Instructor Format
    %
    % Description:
    % The Dhole Optimization Algorithm (DOA) is a nature-inspired metaheuristic
    % based on the cooperative hunting behavior of dholes. The algorithm
    % consists of searching, encircling, and hunting stages to balance
    % exploration and exploitation during optimization.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of search agents)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Best solution vector found
    %   convergenceCurve - Best fitness value at each iteration
    %
    % Tunable Parameters (internal):
    %   C        - Control parameter decreasing linearly (exploration factor)
    %   PWN      - Prey weight number (random integer in [5,20])
    % -------------------------------------------------------------------------
    % Structure:
    %   1. Initialization
    %   2. Fitness Evaluation
    %   3. Main Optimization Loop
    %   4. Boundary Control
    %   5. Best Solution Update
    %   6. Convergence Curve Storage
    % -------------------------------------------------------------------------

    %% --------------------------- Initialization ----------------------------
    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end

    X = zeros(nPop, dim);
    for d = 1:dim
        X(:, d) = rand(nPop, 1) .* (ub(d) - lb(d)) + lb(d);
    end

    fitness = zeros(nPop, 1);
    bestFitness = inf;

    %% ------------------------ Fitness Evaluation ---------------------------
    for i = 1:nPop
        fitness(i) = objFun(X(i, :));
        if fitness(i) < bestFitness
            bestFitness = fitness(i);
            bestPosition = X(i, :);
        end
    end

    localBest_position = bestPosition;
    convergenceCurve = zeros(1, maxItr);

    %% --------------------------- Main Loop ---------------------------------
    for t = 1:maxItr

        C = 1 - (t / maxItr);                 % Eq.(7)
        PWN = round(rand * 15 + 5);           % Eq.(3)
        prey = (bestPosition + localBest_position) / 2;  % Eq.(5)
        prey_local = localBest_position;

        Xnew = zeros(nPop, dim);

        for i = 1:nPop

            if rand < 0.5
                if PWN < 10
                    % ------------------ Searching Stage ---------------------
                    Xnew(i, :) = X(i, :) + C * rand(1, dim) .* (prey - X(i, :));  % Eq.(6)
                else
                    % ------------------ Encircling Stage --------------------
                    for j = 1:dim
                        z = randi([1 nPop]);
                        while z == i
                            z = randi([1 nPop]);
                        end
                        Xnew(i, j) = X(i, j) - X(z, j) + prey(j);  % Eq.(8)
                    end
                end
            else
                % ---------------------- Hunting Stage ------------------------
                Q = 3 * rand * fitness(i) / objFun(prey_local);   % Eq.(10)

                if Q > 2
                    W_prey = exp(-1 / Q) .* prey_local;            % Eq.(11)
                    for j = 1:dim
                        factor = p_obj(PWN);
                        Xnew(i, j) = X(i, j) + ...
                            cos(2 * pi * rand) * W_prey(j) * factor - ...
                            sin(2 * pi * rand) * W_prey(j) * factor;  % Eq.(12)
                    end
                else
                    factor = p_obj(PWN);
                    Xnew(i, :) = (X(i, :) - bestPosition) * factor + ...
                        factor .* rand(1, dim) .* X(i, :);  % Eq.(13)
                end
            end
        end

        %% ---------------------- Boundary Enforcement -----------------------
        for i = 1:nPop
            Xnew(i, :) = min(Xnew(i, :), ub);
            Xnew(i, :) = max(Xnew(i, :), lb);
        end

        %% ------------------ Local & Global Best Update ---------------------
        localBest_position = Xnew(1, :);
        localBest_fitness = objFun(localBest_position);

        for i = 1:nPop
            newFitness = objFun(Xnew(i, :));

            % Local best
            if newFitness < localBest_fitness
                localBest_fitness = newFitness;
                localBest_position = Xnew(i, :);
            end

            % Population update
            if newFitness < fitness(i)
                fitness(i) = newFitness;
                X(i, :) = Xnew(i, :);

                % Global best
                if fitness(i) < bestFitness
                    bestFitness = fitness(i);
                    bestPosition = X(i, :);
                end
            end
        end

        %% -------------------- Convergence Curve ----------------------------
        convergenceCurve(t) = bestFitness;

    end

end

%% ========================== Helper Function ============================
function y = p_obj(x)
    % Eq.(4) Probability-related function
    y = ((1 / (1 + exp(-0.5 * (x - 25))))^2) * rand;
end
