function [Best_FF, Best_P, conv] = RLQFAO(lb, ub, dim, pop, Max_Iter, fobj)
    %% RLQFAO: Reinforced Levy Quantum Fish Algorithm Optimizer
    %  Author: Original RLQFAO code by unknown
    %
    %  Description:
    %    RLQFAO is a metaheuristic optimizer combining Levy flight,
    %    quantum-inspired moves, and reinforcement learning to find
    %    the global minimum of an objective function.
    %
    %  Inputs:
    %    pop       : Number of agents (integer)
    %    Max_Iter  : Maximum number of iterations (integer)
    %    lb        : Lower bounds (1 x dim vector or scalar)
    %    ub        : Upper bounds (1 x dim vector or scalar)
    %    dim       : Problem dimension (integer)
    %    fobj      : Objective function handle
    %
    %  Outputs:
    %    Best_FF   : Best fitness value found
    %    Best_P    : Best solution position (1 x dim)
    %    conv      : Convergence curve (1 x Max_Iter)
    %
    %  Tunable parameters:
    %    alpha, delta : small constants for initialization
    %    Q_table      : reinforcement learning table (state-action)
    %    gamma        : discount factor for RL

    %% --------------------- Initialization ---------------------
    Best_P  = zeros(1, dim);
    Best_FF = inf;

    % Initialize population using Goodnode strategy
    X     = Goodnode_initial(pop, ub, lb, dim);
    Xnew  = X;

    Ffun       = inf(1, pop);
    Ffun_new   = inf(1, pop);

    for i = 1:pop
        Ffun(i) = fobj(X(i,:));
        if Ffun(i) < Best_FF
            Best_FF = Ffun(i);
            Best_P  = X(i,:);
        end
    end

    t      = 1;
    alpha  = 0.1;
    delta  = 0.1;

    action_num = 3;
    state_num  = 3;
    cur_state  = 2;   % initial state
    episode    = 300;

    Reward_table = [1 1 1; 1 1 -1; 1 1 1];
    Q_table      = zeros(state_num, action_num);
    gamma        = 0.8;

    conv = zeros(1, Max_Iter);

    %% --------------------- Main Optimization Loop ---------------------
    while t <= Max_Iter

        % Time-dependent parameters
        G2 = 2*rand()-1;
        G1 = 2*(1 - t/Max_Iter);
        to = 1:dim;
        u  = 0.0265;
        r0 = 10;
        r  = r0 + u*to;
        omega  = 0.005;
        phi0   = 3*pi/2;
        phi    = -omega*to + phi0;
        x_coord = r .* sin(phi);
        y_coord = r .* cos(phi);
        QF = t^((2*rand()-1)/(1-Max_Iter)^2);

        % Strategy selection
        riMOA = 0.5 - 0.5*sin(t/Max_Iter*pi - pi/2);

        for i = 1:pop
            %% --------------------- Strategy 1 & 2 ---------------------
            if rand <= riMOA
                if rand < 0.5
                    % Quantum-inspired rotational moves
                    alpha1 = zeros(1,dim);
                    beta   = zeros(1,dim);
                    for j = 1:dim
                        alpha1(j) = Best_P(j)/norm(Best_P);
                        if rand < 0.5
                            P = ones(1,dim);
                        else
                            P = -ones(1,dim);
                        end
                        beta(j) = P(j) * sqrt(1 - (Best_P(j)/norm(Best_P))^2);
                    end

                    X_tema = alpha1 .* norm(Best_P);
                    X_temb = beta   .* norm(Best_P);
                    Z_R = [alpha1; beta];
                    theta = 2*pi*rand;
                    R_theta = [cos(theta), -sin(theta); sin(theta), cos(theta)];
                    Z_R2 = R_theta * Z_R;
                    X_rtema = Z_R2(1,:) .* norm(Best_P);
                    X_rtemb = Z_R2(2,:) .* norm(Best_P);

                    x_b1 = Best_P.*X_tema + (mean(X(i,:))-Best_P).*randn(1,dim);
                    x_b2 = Best_P.*X_temb + (mean(X(i,:))-Best_P).*randn(1,dim);
                    x_b3 = Best_P.*X_rtema + (mean(X(i,:))-Best_P).*randn(1,dim);
                    x_b4 = Best_P.*X_rtemb + (mean(X(i,:))-Best_P).*randn(1,dim);

                    x_b1 = Bounds(x_b1, lb, ub);
                    x_b2 = Bounds(x_b2, lb, ub);
                    x_b3 = Bounds(x_b3, lb, ub);
                    x_b4 = Bounds(x_b4, lb, ub);

                    fitness = [fobj(x_b1), fobj(x_b2), fobj(x_b3), fobj(x_b4)];
                    x_b = [x_b1; x_b2; x_b3; x_b4];
                    [F_X2, index] = min(fitness);
                    Xnew(i,:) = x_b(index,:);

                    if F_X2 < Ffun(i)
                        X(i,:) = Xnew(i,:);
                        Ffun(i) = F_X2;
                    end

                else
                    % Levy flight & spiral move
                    Xnew(i,:) = Best_P.*Levy(dim) + X(floor(pop*rand+1),:) + (y_coord - x_coord)*rand;
                    Xnew(i,:) = Bounds(Xnew(i,:), lb, ub);
                    F_X2 = fobj(Xnew(i,:));
                    if F_X2 < Ffun(i)
                        X(i,:) = Xnew(i,:);
                        Ffun(i) = F_X2;
                    end
                end

            else
                %% --------------------- Strategy 3 ---------------------
                Flame_no = round(pop - t*((pop-1)/Max_Iter));
                a = -1 + t*((-1)/Max_Iter);
                [bb, aa] = sort(Ffun);

                if rand < 0.5
                    for j = 1:dim
                        if i <= Flame_no
                            distance = abs(Best_P(j)-X(i,j));
                            b = 1.5;
                            Tt = (a-1)*rand +1;
                            Xnew(i,j) = distance*exp(b*Tt)*cos(Tt*2*pi) + Best_P(j);
                        else
                            distance = abs(Best_P(j)-X(i,j));
                            b = 1.5;
                            Tt = (a-1)*rand +1;
                            Xnew(i,j) = distance*exp(b*Tt)*cos(Tt*2*pi) + X(aa(Flame_no),j);
                        end
                    end
                else
                    Xnew(i,:) = QF*Best_P - (G2*X(i,:)*rand) - G1*Levy(dim) + rand*G2;
                end

                Xnew(i,:) = Bounds(Xnew(i,:), lb, ub);
                F_X2 = fobj(Xnew(i,:));
                if F_X2 < Ffun(i)
                    X(i,:) = Xnew(i,:);
                    Ffun(i) = F_X2;
                end
            end

            %% --------------------- Strategy 4 & Reinforcement Learning ---------------------
            if (Q_table(cur_state,1) >= Q_table(cur_state,2)) && (Q_table(cur_state,1) >= Q_table(cur_state,3))
                action_num = 1;
                k = 12000;
                if rand < 0.5
                    Xnew(i,:) = lb + ub - rand*X(i,:);
                else
                    Xnew(i,:) = (ub+lb)/2 + (ub+lb)/(2*k) - X(i,:)/k;
                end
            elseif (Q_table(cur_state,2) >= Q_table(cur_state,1)) && (Q_table(cur_state,2) >= Q_table(cur_state,3))
                action_num = 2;
                F = normrnd(0.3,0.1);
                Xnew1 = 0.1*X(i,:) + 1.5*F*(mean(X(:,:)) - X(i,:)) + 1.5*F*(Best_P - X(i,:));
                Xnew(i,:) = 0.8*X(i,:) + Xnew1;
            else
                action_num = 3;
                Xnew(i,:) = X(i,:) + rand*(Best_P - X(i,:)).*Levy(dim);
            end

            Xnew(i,:) = Bounds(Xnew(i,:), lb, ub);
            Ffun_new(i) = fobj(Xnew(i,:));

            if Ffun_new(i) <= Ffun(i)
                Reward_table(cur_state, action_num) = +1;
                X(i,:) = Xnew(i,:);
                Ffun(i) = Ffun_new(i);
            else
                Reward_table(cur_state, action_num) = -1;
            end
            cur_state = action_num;

            %% --------------------- Q-table Update ---------------------
            current_state = randi(state_num,1);
            for k = 1:episode
                while isempty(find(Reward_table(current_state,:) > -1,1))
                    current_state = rem(current_state, state_num)+1;
                end
                possible_actions = find(Reward_table(current_state,:) > -1);
                chosen_action = possible_actions(randi(length(possible_actions)));
                r = Reward_table(current_state, chosen_action);
                next_state = chosen_action;
                if isempty(find(Reward_table(next_state,:) > -1,1))
                    next_state = rem(next_state, state_num)+1;
                end
                next_actions = find(Reward_table(next_state,:) > -1);
                maxQ = max(Q_table(next_state, next_actions));
                study_rate = 1 - 0.9*t/Max_Iter;
                Q_table(current_state, chosen_action) = Q_table(current_state, chosen_action) + ...
                    study_rate*(r + gamma*maxQ - Q_table(current_state, chosen_action));
                current_state = next_state;
            end
        end

        %% --------------------- Best Solution Update ---------------------
        for i = 1:pop
            if Ffun(i) < Best_FF
                Best_FF = Ffun(i);
                Best_P  = X(i,:);
            end
        end

        conv(t) = Best_FF;
        t = t+1;
    end

end

%% --------------------- Helper Functions ---------------------
function o = Levy(d)
    beta = 1.5;
    sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u = randn(1,d)*sigma; v = randn(1,d);
    o = u ./ abs(v).^(1/beta);
end

function s = Bounds(s, Lb, Ub)
    a = rand(1,length(Lb));
    s(s < Lb) = Lb(s < Lb) + (Ub(s < Lb) - Lb(s < Lb)).*a(s < Lb);
    s(s > Ub) = Lb(s > Ub) + (Ub(s > Ub) - Lb(s > Ub)).*a(s > Ub);
end

function X = Goodnode_initial(pop, ub, lb, dim)
    GD = Goodnode(pop, dim);
    X = zeros(pop, dim);
    for i = 1:pop
        for j = 1:dim
            X(i,j) = (ub(j)-lb(j))*GD(i,j) + lb(j);
        end
    end
end

function GD = Goodnode(M,N)
    if nargin == 0
        M = 1000; N = 2;
    end
    tmp1 = [1:M]'*ones(1,N);
    Ind = 1:N;
    prime1 = primes(100*N);
    [p,q] = find(prime1 >= (2*N+3));
    tmp2 = (2*pi.*Ind)/prime1(1,q(1));
    tmp2 = 2*cos(tmp2);
    tmp2 = ones(M,1)*tmp2;
    GD = mod(tmp1.*tmp2,1);
end
