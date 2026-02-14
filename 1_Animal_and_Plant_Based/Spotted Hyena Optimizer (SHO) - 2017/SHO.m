function [Best_hyena_score, Best_hyena_pos, Convergence_curve] = SHO(lb, ub, dim, nPop, maxItr, objFun)
    %__________________________________________________________________________
    % Social Hyena Optimizer (SHO)
    %__________________________________________________________________________
    % Author: Original code by user (refactored)
    %
    % Description:
    %   This function implements the Social Hyena Optimizer (SHO) for
    %   unconstrained minimization problems.
    %
    % Inputs:
    %   lb       - Lower bound (scalar or 1xdim vector)
    %   ub       - Upper bound (scalar or 1xdim vector)
    %   dim      - Problem dimensionality
    %   nPop     - Number of hyena agents (population size)
    %   maxItr   - Maximum number of iterations
    %   objFun   - Handle to the objective function: fitness = objFun(position)
    %
    % Outputs:
    %   Best_hyena_score - Best fitness value found
    %   Best_hyena_pos   - Position of the best solution
    %   Convergence_curve- Best fitness at each iteration (1 x maxItr)
    %
    % Tunable Parameters:
    %   a - Exploration-exploitation control, decreases linearly from 5 to 0

    %% Initialization
    hyena_pos = initPopulation(nPop, dim, ub, lb);
    Convergence_curve = zeros(1, maxItr);

    pre_population = hyena_pos;
    pre_fitness = zeros(1, nPop);

    Iteration = 1;

    %% Main Loop
    while Iteration <= maxItr
        % Boundary check and fitness evaluation
        hyena_fitness = zeros(1, nPop);
        for i = 1:nPop
            hyena_pos(i,:) = min(max(hyena_pos(i,:), lb), ub);
            hyena_fitness(i) = objFun(hyena_pos(i,:));
        end

        % Update best hyenas
        if Iteration == 1
            [fitness_sorted, FS] = sort(hyena_fitness);
            sorted_population = hyena_pos(FS,:);
            best_hyenas = sorted_population;
            best_hyena_fitness = fitness_sorted;
        else
            double_population = [pre_population; best_hyenas];
            double_fitness = [pre_fitness, best_hyena_fitness];
            [double_fitness_sorted, FS] = sort(double_fitness);
            double_sorted_population = double_population(FS,:);
            fitness_sorted = double_fitness_sorted(1:nPop);
            sorted_population = double_sorted_population(1:nPop,:);
            best_hyenas = sorted_population;
            best_hyena_fitness = fitness_sorted;
        end

        % Compute number of influencing hyenas
        NOH = computeNOH(best_hyena_fitness);

        % Update global best
        Best_hyena_score = fitness_sorted(1);
        Best_hyena_pos = sorted_population(1,:);

        pre_population = hyena_pos;
        pre_fitness = hyena_fitness;

        % Adaptive control parameter
        a = 5 - Iteration * (5 / maxItr);

        % Position update
        for i = 1:nPop
            for j = 1:dim
                CV = 0;
                for k = 1:NOH
                    r1 = rand();
                    r2 = rand();
                    Var1 = 2 * a * r1 - a;
                    Var2 = 2 * r2;
                    distance_to_hyena = abs(Var2 * sorted_population(k,j) - hyena_pos(i,j));
                    CV = CV + (sorted_population(k,j) - Var1 * distance_to_hyena);
                end
                hyena_pos(i,j) = CV / (NOH + 1);
            end
        end

        % Store convergence
        Convergence_curve(Iteration) = Best_hyena_score;
        Iteration = Iteration + 1;
    end

end

%% ------------------ Helper Functions ------------------ %%
function X = computeNOH(best_hyena_fitness)
    % Determine number of influencing hyenas
    min_val = 0.5; max_val = 1;
    M = (max_val - min_val) * rand() + min_val + best_hyena_fitness(1);
    count = sum(M >= best_hyena_fitness(2:end));
    X = count;
end

function Pos = initPopulation(SearchAgents, dim, ub, lb)
    % Initialize population within bounds
    if isscalar(ub)
        Pos = rand(SearchAgents, dim) .* (ub - lb) + lb;
    else
        Pos = zeros(SearchAgents, dim);
        for i = 1:dim
            Pos(:,i) = rand(SearchAgents,1) .* (ub(i) - lb(i)) + lb(i);
        end
    end
end
