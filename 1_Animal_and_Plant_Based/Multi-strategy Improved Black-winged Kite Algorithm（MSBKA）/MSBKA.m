function [Best_Fitness_MSBKA, Best_Pos_MSBKA, Convergence_curve] = MSBKA(lb, ub, dim, pop, T, fobj)
    %% ------------------------------------------------------------------------
    % Multi-strategy Improved Black-winged Kite Algorithm (MSBKA)
    % -------------------------------------------------------------------------
    % Author: Original algorithm source (adapted for MATLAB benchmark)
    %
    % Description:
    %   MSBKA is a population-based metaheuristic for global optimization.
    %   It combines Levy-flight global search, differential mutation-based
    %   local exploitation, and migration strategies guided by the current best.
    %
    % Inputs:
    %   T       - Maximum number of iterations (integer)
    %   pop     - Population size (integer)
    %   lb      - Lower bounds (scalar or 1xdim vector)
    %   ub      - Upper bounds (scalar or 1xdim vector)
    %   dim     - Problem dimensionality (integer)
    %   fobj    - Objective function handle (fitness to minimize)
    %
    % Outputs:
    %   Best_Fitness_MSBKA - Best fitness value found
    %   Best_Pos_MSBKA     - Position vector of best solution
    %   Convergence_curve  - 1xT vector of best fitness per iteration
    %
    % Tunable parameters:
    %   p       - Attack probability (linearly decreases from 0.9 to 0)
    %   n       - Step scaling factor for global search (dynamic)
    %   F       - Differential mutation factor (dynamic)
    % -------------------------------------------------------------------------

    %% ------------------ Initialization -------------------------------------
    % Ensure bounds are vectors
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    % Tent chaotic initialization
    XPos = Tent_initialization(pop, dim, lb, ub);

    % Evaluate initial fitness
    XFit = zeros(pop,1);
    for i = 1:pop
        XFit(i) = fobj(XPos(i,:));
    end

    Convergence_curve = zeros(1, T);
    Best_Fitness_MSBKA = inf;
    Best_Pos_MSBKA = zeros(1, dim);

    %% ------------------ Main Iteration Loop --------------------------------
    for t = 1:T
        % Sort by fitness
        [XFit_sorted, sorted_idx] = sort(XFit);
        XLeader_Pos = XPos(sorted_idx(1), :);
        XLeader_Fit = XFit_sorted(1);

        % Dynamic attack probability
        p = 0.9 * (1 - t/T);

        XPosNew = XPos; % Pre-allocate
        XFit_New = zeros(pop,1);

        for i = 1:pop
            n = 0.05 * exp(-2 * (t/T)^2);
            levy_step = levyFlight(dim);
            cauchy_step = trnd(1, 1, dim);

            if rand > p
                % Strategy 1: Levy-flight guided global search
                XPosNew(i,:) = XPos(i,:) + n .* levy_step .* (XLeader_Pos - XPos(i,:));
            else
                % Strategy 2: Differential mutation + local exploitation
                r_idx = randperm(pop, 3);
                F = 0.5 + 0.3 * sin(pi*t/T);
                XPosNew(i,:) = XPos(r_idx(1),:) + F * (XPos(r_idx(2),:) - XPos(r_idx(3),:));
                XPosNew(i,:) = XPosNew(i,:) + 0.1 * cauchy_step .* (XLeader_Pos - XPosNew(i,:));
            end

            % Reflective boundary handling
            XPosNew(i,:) = reflectBounds(XPosNew(i,:), lb, ub);

            % Evaluate fitness
            XFit_New(i) = fobj(XPosNew(i,:));

            % Greedy selection
            if XFit_New(i) < XFit(i)
                XPos(i,:) = XPosNew(i,:);
                XFit(i) = XFit_New(i);
            end

            %% --------------- Migration Behavior --------------------
            m = 2 * sin(rand + pi/2);
            golden_ratio = 0.5*(sqrt(5)-1);
            x1 = -pi + (1 - golden_ratio) * 2*pi;
            x2 = -pi + golden_ratio * 2*pi;

            R1 = rand() * pi;
            R2 = rand() * 2*pi;
            s = randi([1, pop]);
            r_XFitness = XFit(s);
            ori_value = rand(1, dim);
            cauchy_value = tan((ori_value - 0.5) * pi);

            if XFit(i) < r_XFitness
                XPosNew(i,:) = XPos(i,:) + R2 * sin(R1) .* cauchy_value .* (x1*XPos(i,:) - x2*XLeader_Pos);
            else
                XPosNew(i,:) = XPos(i,:) + cauchy_value .* (XLeader_Pos - m .* XPos(i,:));
            end

            % Boundary enforcement
            XPosNew(i,:) = max(min(XPosNew(i,:), ub), lb);

            % Evaluate fitness after migration
            XFit_New(i) = fobj(XPosNew(i,:));

            % Greedy selection
            if XFit_New(i) < XFit(i)
                XPos(i,:) = XPosNew(i,:);
                XFit(i) = XFit_New(i);
            end
        end

        %% ----------------- Update Global Best ----------------------
        [current_best_fit, idx] = min(XFit);
        if current_best_fit < Best_Fitness_MSBKA
            Best_Fitness_MSBKA = current_best_fit;
            Best_Pos_MSBKA = XPos(idx,:);
        end

        Convergence_curve(t) = Best_Fitness_MSBKA;
    end

end

%% =================== Auxiliary Functions ===============================

function positions = reflectBounds(positions, lb, ub)
    % Reflective boundary handling
    for j = 1:length(positions)
        if positions(j) < lb(j)
            positions(j) = 2 * lb(j) - positions(j);
        elseif positions(j) > ub(j)
            positions(j) = 2 * ub(j) - positions(j);
        end
        positions(j) = min(max(positions(j), lb(j)), ub(j));
    end
end

function step = levyFlight(d)
    % Levy flight step generation
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
        (gamma((1 + beta)/2) * beta * 2^((beta - 1)/2)))^(1/beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1/beta);
end

function Positions = Tent_initialization(SearchAgents_no, dim, lb, ub)
    % Tent chaotic map initialization
    Positions = zeros(SearchAgents_no, dim);
    epsilon = 1e-10;
    for i = 1:SearchAgents_no
        x = rand;
        while abs(x - 0.5) < epsilon, x = rand; end
        for j = 1:dim
            if x < 0.5
                x = 2 * x;
            else
                x = 2 * (1 - x);
            end
            x = min(max(x, epsilon), 1 - epsilon);
            Positions(i,j) = lb(j) + x * (ub(j) - lb(j));
        end
    end
end
