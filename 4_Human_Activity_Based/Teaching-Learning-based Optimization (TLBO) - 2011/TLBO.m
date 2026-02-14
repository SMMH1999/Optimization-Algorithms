function [bestFitness, bestPosition, convergenceCurve] = TLBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Teaching-Learning-Based Optimization (TLBO)
    % -------------------------------------------------------------------------
    % Developer   : S. Mostapha Kalami Heris (Yarpiz Team)
    % Refactored  : Benchmark unified format (Instructor Style)
    %
    % Description:
    % TLBO is a population-based metaheuristic inspired by the teaching and
    % learning process in a classroom. The algorithm consists of two main
    % phases: Teacher Phase (global guidance) and Learner Phase (peer learning).
    % The best individual acts as the teacher and guides the population toward
    % better solutions. The algorithm assumes minimization.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (number of learners)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (fitness = objFun(position))
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (maxItr×1)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   nPop   - Controls population diversity
    %   maxItr - Controls termination criterion
    %
    % Notes:
    %   • Supports scalar or vector bounds.
    %   • Boundary constraints are strictly enforced.
    %   • No function evaluation counters (FEs) included.
    %   • Designed for benchmark optimization frameworks.
    % =========================================================================

    %% --------------------------- Initialization -----------------------------

    % Ensure bounds are row vectors
    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    % Initialize population
    pop = zeros(nPop, dim);
    fitness = zeros(nPop, 1);

    for i = 1:nPop
        pop(i, :) = lb + rand(1, dim) .* (ub - lb);
        fitness(i) = objFun(pop(i, :));
    end

    % Global best initialization
    [bestFitness, idx] = min(fitness);
    bestPosition = pop(idx, :);

    % Convergence history
    convergenceCurve = zeros(maxItr, 1);

    %% ============================ Main Loop =================================
    for t = 1:maxItr

        %% ------------------------- Teacher Phase ----------------------------

        % Compute mean of population
        Mean = mean(pop, 1);

        % Identify teacher (best individual)
        [~, teacherIdx] = min(fitness);
        Teacher = pop(teacherIdx, :);

        for i = 1:nPop

            % Teaching Factor (either 1 or 2)
            TF = randi([1 2]);

            % Generate new solution
            newPosition = pop(i, :) + rand(1, dim) .* (Teacher - TF * Mean);

            % Boundary enforcement
            newPosition = max(newPosition, lb);
            newPosition = min(newPosition, ub);

            % Evaluate
            newFitness = objFun(newPosition);

            % Greedy selection
            if newFitness < fitness(i)
                pop(i, :) = newPosition;
                fitness(i) = newFitness;

                % Update global best
                if newFitness < bestFitness
                    bestFitness = newFitness;
                    bestPosition = newPosition;
                end
            end
        end

        %% ------------------------- Learner Phase ----------------------------

        for i = 1:nPop

            % Select random partner j ≠ i
            j = randi(nPop);
            while j == i
                j = randi(nPop);
            end

            Step = pop(i, :) - pop(j, :);
            if fitness(j) < fitness(i)
                Step = -Step;
            end

            % Generate new solution
            newPosition = pop(i, :) + rand(1, dim) .* Step;

            % Boundary enforcement
            newPosition = max(newPosition, lb);
            newPosition = min(newPosition, ub);

            % Evaluate
            newFitness = objFun(newPosition);

            % Greedy selection
            if newFitness < fitness(i)
                pop(i, :) = newPosition;
                fitness(i) = newFitness;

                % Update global best
                if newFitness < bestFitness
                    bestFitness = newFitness;
                    bestPosition = newPosition;
                end
            end
        end

        %% ---------------------- Convergence Update --------------------------

        convergenceCurve(t) = bestFitness;

    end

end
