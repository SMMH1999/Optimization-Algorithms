function [Best_Fitness, Best_Position, Convergence_curve] = WOA(LB, UB, Dim, SearchAgents_no, Max_iter, Cost_Function)
    % WOA - Whale Optimization Algorithm
    % Algorithm constant
    b = 1;

    % Expand bounds if scalar
    if isscalar(UB)
        UB = repmat(UB, 1, Dim);
        LB = repmat(LB, 1, Dim);
    end

    % Initialization
    Positions = initialization(SearchAgents_no, Dim, UB, LB);
    Best_Fitness  = inf;
    Best_Position = zeros(1, Dim);
    Convergence_curve = inf(1, Max_iter);

    % Main loop
    for t = 1:Max_iter

        % Linearly decreasing parameter a
        a = 2 - t * (2 / Max_iter);

        %% Fitness evaluation and leader update
        for i = 1:SearchAgents_no

            % Boundary control
            Positions(i,:) = max(min(Positions(i,:), UB), LB);

            % Fitness
            fit = Cost_Function(Positions(i,:));

            if fit < Best_Fitness
                Best_Fitness  = fit;
                Best_Position = Positions(i,:);
            end
        end

        Convergence_curve(t) = Best_Fitness;

        %% Update whale positions
        for i = 1:SearchAgents_no

            r1 = rand();
            r2 = rand();
            A  = 2 * a * r1 - a;
            C  = 2 * r2;
            l  = -1 + 2 * rand();   % in [-1,1]
            p  = rand();

            for j = 1:Dim

                if p < 0.5
                    if abs(A) >= 1
                        % Exploration: random whale
                        rand_idx = randi(SearchAgents_no);
                        X_rand   = Positions(rand_idx, j);
                        D_rand   = abs(C * X_rand - Positions(i, j));
                        Positions(i, j) = X_rand - A * D_rand;
                    else
                        % Exploitation: encircle best
                        D_best = abs(C * Best_Position(j) - Positions(i, j));
                        Positions(i, j) = Best_Position(j) - A * D_best;
                    end
                else
                    % Exploitation: spiral motion
                    D_best = abs(Best_Position(j) - Positions(i, j));
                    Positions(i, j) = D_best * exp(b * l) * cos(2*pi*l) ...
                        + Best_Position(j);
                end
            end
        end
    end
end

%% ====================================================================

% Population initialization
function X = initialization(N, dim, UB, LB)
    X = LB + rand(N, dim) .* (UB - LB);
end
