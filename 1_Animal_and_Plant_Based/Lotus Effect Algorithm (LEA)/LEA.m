function [Best_score, Best_pos, cg_curve] = LEA(lb, ub, dim, SearchAgents_no, Max_iteration, Cost_Function)

    cg_curve = zeros(1, Max_iteration);

    % Expand bounds if scalar
    if isscalar(ub)
        ub = ones(1, dim) * ub;
        lb = ones(1, dim) * lb;
    end

    Delta_max = (ub - lb) / 10;

    % Best (food) and worst (enemy)
    Food_fitness  = inf;
    Food_pos      = zeros(1, dim);
    Enemy_fitness = -inf;
    Enemy_pos     = zeros(1, dim);

    % Initialize population
    X      = Population_Generator(SearchAgents_no, dim, ub, lb);
    DeltaX = Population_Generator(SearchAgents_no, dim, ub, lb);
    Fitness = zeros(1, SearchAgents_no);

    for iter = 1:Max_iteration

        % Dynamic parameters
        r = (ub - lb)/4 + (ub - lb) * (iter/Max_iteration) * 2;
        w = 0.9 - iter * ((0.9 - 0.4) / Max_iteration);

        my_c = 0.1 - iter * (0.1 / (Max_iteration/2));
        my_c = max(my_c, 0);

        s = 2 * rand * my_c;
        a = 2 * rand * my_c;
        c = 2 * rand * my_c;
        f = 2 * rand;
        e = my_c;

        %% Fitness evaluation
        for i = 1:SearchAgents_no
            Fitness(i) = Cost_Function(X(i, :));

            % Update food
            if Fitness(i) < Food_fitness
                Food_fitness = Fitness(i);
                Food_pos = X(i, :);
            end

            % Update enemy
            if Fitness(i) > Enemy_fitness
                if all(X(i, :) < ub) && all(X(i, :) > lb)
                    Enemy_fitness = Fitness(i);
                    Enemy_pos = X(i, :);
                end
            end
        end

        %% Update agents
        for i = 1:SearchAgents_no

            neighbours_no = 0;
            Neighbours_X = [];
            Neighbours_DeltaX = [];

            % Find neighbors
            for j = 1:SearchAgents_no
                Dist = distance(X(i, :), X(j, :));
                if all(Dist <= r) && any(Dist ~= 0)
                    neighbours_no = neighbours_no + 1;
                    Neighbours_X(neighbours_no, :) = X(j, :);
                    Neighbours_DeltaX(neighbours_no, :) = DeltaX(j, :);
                end
            end

            % Separation
            if neighbours_no > 1
                S = -sum(Neighbours_X - X(i, :), 1);
                A = mean(Neighbours_DeltaX, 1);
                C_temp = mean(Neighbours_X, 1);
            else
                S = zeros(1, dim);
                A = DeltaX(i, :);
                C_temp = X(i, :);
            end

            % Cohesion
            C = C_temp - X(i, :);

            % Attraction to food
            Dist2Food = distance(X(i, :), Food_pos);
            if all(Dist2Food <= r)
                F = Food_pos - X(i, :);
            else
                F = zeros(1, dim);
            end

            % Distraction from enemy
            Dist2Enemy = distance(X(i, :), Enemy_pos);
            if all(Dist2Enemy <= r)
                Enemy = Enemy_pos + X(i, :);
            else
                Enemy = zeros(1, dim);
            end

            % Boundary handling
            for d = 1:dim
                if X(i, d) > ub(d)
                    X(i, d) = lb(d);
                    DeltaX(i, d) = rand;
                elseif X(i, d) < lb(d)
                    X(i, d) = ub(d);
                    DeltaX(i, d) = rand;
                end
            end

            %% Position update
            if any(Dist2Food > r)
                if neighbours_no > 1
                    for d = 1:dim
                        DeltaX(i, d) = w * DeltaX(i, d) + rand * A(d) + rand * C(d) + rand * S(d);
                        DeltaX(i, d) = max(min(DeltaX(i, d), Delta_max(d)), -Delta_max(d));
                        X(i, d) = X(i, d) + DeltaX(i, d);
                    end
                else
                    % Levy flight
                    X(i, :) = X(i, :) + Levy(dim) .* X(i, :);
                    DeltaX(i, :) = 0;
                end
            else
                for d = 1:dim
                    DeltaX(i, d) = (a*A(d) + c*C(d) + s*S(d) + f*F(d) + e*Enemy(d)) + w * DeltaX(i, d);
                    DeltaX(i, d) = max(min(DeltaX(i, d), Delta_max(d)), -Delta_max(d));
                    X(i, d) = X(i, d) + DeltaX(i, d);
                end
            end

            % Final bounds check
            Flag4ub = X(i, :) > ub;
            Flag4lb = X(i, :) < lb;
            X(i, :) = X(i, :) .* ~(Flag4ub + Flag4lb) + ub .* Flag4ub + lb .* Flag4lb;
        end

        % Store best
        Best_score = Food_fitness;
        Best_pos   = Food_pos;
        cg_curve(iter) = Best_score;
    end
end

%% ------------------ Helper Functions ------------------- %%
function o = distance(a,b)
    for i=1:size(a,1)
        o(1,i)=sqrt((a(i)-b(i))^2);
    end
end

function o = Euclidean_distance(a,b)
    for i=1:size(a,1)
        o(1,i)=sqrt((a(i)-b(i))^2);
    end
end

function o=Levy(d)
    beta=3/2;

    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;
    v=randn(1,d);
    step=u./abs(v).^(1/beta);

    o=0.01*step;
end