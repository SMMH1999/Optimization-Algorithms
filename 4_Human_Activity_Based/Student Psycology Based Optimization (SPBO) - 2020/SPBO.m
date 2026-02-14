function [bestFitness, bestPosition, convergenceCurve] = SPBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Student Psychology Based Optimization (SPBO)
    % Abbreviation: SPBO
    %
    % Reference:
    % Bikash Das, V. Mukherjee, D. Das,
    % "Student psychology based optimization algorithm: A new population based
    % optimization algorithm for solving optimization problems",
    % Advances in Engineering Software, 146 (2020) 102804.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (1×dim vector or scalar)
    %   ub        - Upper bound (1×dim vector or scalar)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Number of students (population size)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Best solution vector (1×dim)
    %   convergenceCurve - Best fitness at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   nPop     - Controls population diversity and convergence speed
    %   maxItr   - Controls total search duration
    %
    % Notes:
    %   - Minimization problem assumed.
    %   - Fully benchmark-compatible format.
    %   - No display or external dependencies.
    % =========================================================================

    %% --------------------------- Initialization -----------------------------

    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    % Initialize population
    X = rand(nPop, dim) .* (ub - lb) + lb;

    sol = X;                          % Current accepted solutions
    ansFit = zeros(nPop, 1);          % Fitness values

    % Evaluate initial population
    for i = 1:nPop
        ansFit(i) = objFun(X(i, :));
    end

    % Determine global best
    [bestFitness, idx] = min(ansFit);
    bestPosition = sol(idx, :);

    convergenceCurve = zeros(1, maxItr);

    %% ------------------------------ Main Loop -------------------------------

    for t = 1:maxItr

        for d = 1:dim

            % Compute mean of each dimension
            meanSol = mean(sol, 1);

            par = sol;
            par1 = sol;

            check = rand(nPop,1);
            mid   = rand(nPop,1);

            for i = 1:nPop

                if ansFit(i) == bestFitness
                    % Best student (Equation 1)
                    randIdx = randi(nPop);
                    lk = randIdx;
                    par1(i,d) = par(i,d) + ((-1)^(round(1+rand))) * rand * ...
                        (par(i,d) - par(lk,d));

                elseif check(i) < mid(i)
                    % Good student (Equation 2a & 2b)
                    if rand > rand
                        par1(i,d) = bestPosition(d) + ...
                            rand * (bestPosition(d) - par(i,d));
                    else
                        par1(i,d) = par(i,d) + ...
                            rand * (bestPosition(d) - par(i,d)) + ...
                            rand * (par(i,d) - meanSol(d));
                    end
                else
                    % Average / random improving student (Equation 3 & 4)
                    if rand > rand
                        par1(i,d) = par(i,d) + ...
                            rand * (meanSol(d) - par(i,d));
                    else
                        par1(i,d) = lb(d) + rand * (ub(d) - lb(d));
                    end
                end

                % Boundary control
                par1(i,d) = max(par1(i,d), lb(d));
                par1(i,d) = min(par1(i,d), ub(d));

            end

            % Evaluate updated population
            newFit = zeros(nPop,1);
            for i = 1:nPop
                newFit(i) = objFun(par1(i,:));
            end

            % Greedy selection
            for i = 1:nPop
                if ansFit(i) > newFit(i)
                    ansFit(i) = newFit(i);
                    sol(i,:) = par1(i,:);
                end
            end

            % Update global best
            [currentBest, idx] = min(ansFit);
            if currentBest < bestFitness
                bestFitness = currentBest;
                bestPosition = sol(idx,:);
            end

        end

        convergenceCurve(t) = bestFitness;

    end

end
