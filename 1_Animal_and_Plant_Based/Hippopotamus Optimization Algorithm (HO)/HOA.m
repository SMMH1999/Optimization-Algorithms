function [bestFitness, bestPosition, convergenceCurve] = HOA(lb, ub, dim, SearchAgents, Max_iterations, objFun)
    %% Hippopotamus Optimization Algorithm (HOA)
    % Designed and Developed by Mohammad Hussien Amiri and Nastaran Mehrabi Hashjin
    %
    % Description:
    %   HOA is a nature-inspired metaheuristic based on the behavior of
    %   hippopotamuses in rivers and ponds. It combines exploration and
    %   exploitation phases to minimize a given objective function.
    %
    % Input:
    %   SearchAgents   - integer, number of candidate solutions (population size)
    %   Max_iterations - integer, maximum number of iterations
    %   lb             - scalar or 1xD vector, lower bounds of variables
    %   ub             - scalar or 1xD vector, upper bounds of variables
    %   dim            - integer, number of decision variables (dimension)
    %   objFun         - function handle, objective function to minimize
    %
    % Output:
    %   bestFitness       - scalar, best objective function value found
    %   bestPosition      - 1xD vector, position of the best solution
    %   convergenceCurve  - 1xMax_iterations vector, best value per iteration
    %
    % Tunable Parameters:
    %   None beyond the inputs; algorithm uses stochastic parameters internally.

    %% Ensure bounds are vectors
    lb = ones(1, dim) .* lb;
    ub = ones(1, dim) .* ub;

    %% Initialization
    X = lb + rand(SearchAgents, dim) .* (ub - lb);   % population
    fit = zeros(SearchAgents, 1);
    for i = 1:SearchAgents
        fit(i) = objFun(X(i,:));
    end

    best_so_far = zeros(1, Max_iterations);

    %% Main Loop
    for t = 1:Max_iterations
        %% Update the Best Candidate Solution
        [currentBest, idx] = min(fit);
        if t == 1 || currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = X(idx,:);
        end

        %% Phase 1: Exploration (Hippopotamuses in river/pond)
        for i = 1:floor(SearchAgents/2)
            Dominant_hippo = bestPosition;
            I1 = randi([1,2]);
            I2 = randi([1,2]);
            Ip1 = randi([0,1],1,2);

            RandGroupNumber = randi(SearchAgents);
            RandGroup = randperm(SearchAgents, RandGroupNumber);
            MeanGroup = mean(X(RandGroup,:),1) .* (length(RandGroup)~=1) + X(RandGroup(1),:) * (length(RandGroup)==1);

            Alfa = {I2*rand(1,dim)+~Ip1(1), 2*rand(1,dim)-1, rand(1,dim), I1*rand(1,dim)+~Ip1(2), rand(1,dim)};
            A = Alfa{randi([1,5])};
            B = Alfa{randi([1,5])};

            X_P1 = X(i,:) + rand .* (Dominant_hippo - I1 .* X(i,:));

            T = exp(-t/Max_iterations);
            if T > 0.6
                X_P2 = X(i,:) + A .* (Dominant_hippo - I2 .* MeanGroup);
            else
                if rand > 0.5
                    X_P2 = X(i,:) + B .* (MeanGroup - Dominant_hippo);
                else
                    X_P2 = lb + rand(1,dim) .* (ub - lb);
                end
            end

            X_P2 = min(max(X_P2, lb), ub);

            % Evaluate and update
            F_P1 = objFun(X_P1);
            if F_P1 < fit(i)
                X(i,:) = X_P1;
                fit(i) = F_P1;
            end

            F_P2 = objFun(X_P2);
            if F_P2 < fit(i)
                X(i,:) = X_P2;
                fit(i) = F_P2;
            end
        end

        %% Phase 2: Defense against predators (Exploration)
        for i = floor(SearchAgents/2)+1:SearchAgents
            predator = lb + rand(1,dim) .* (ub - lb);
            F_HL = objFun(predator);
            distance2Leader = abs(predator - X(i,:));
            b = unifrnd(2,4);
            c = unifrnd(1,1.5);
            d = unifrnd(2,3);
            l = unifrnd(-2*pi,2*pi);
            RL = 0.05*levy(SearchAgents, dim, 1.5);

            if fit(i) > F_HL
                X_P3 = RL(i,:) .* predator + (b./(c - d*cos(l))) .* (1./distance2Leader);
            else
                X_P3 = RL(i,:) .* predator + (b./(c - d*cos(l))) .* (1./(2*distance2Leader + rand(1,dim)));
            end
            X_P3 = min(max(X_P3, lb), ub);

            F_P3 = objFun(X_P3);
            if F_P3 < fit(i)
                X(i,:) = X_P3;
                fit(i) = F_P3;
            end
        end

        %% Phase 3: Escaping predators (Exploitation)
        for i = 1:SearchAgents
            LO_LOCAL = lb ./ t;
            HI_LOCAL = ub ./ t;
            Alfa = {2*rand(1,dim)-1, rand, randn};
            D = Alfa{randi([1,3])};

            X_P4 = X(i,:) + rand .* (LO_LOCAL + D .* (HI_LOCAL - LO_LOCAL));
            X_P4 = min(max(X_P4, lb), ub);

            F_P4 = objFun(X_P4);
            if F_P4 < fit(i)
                X(i,:) = X_P4;
                fit(i) = F_P4;
            end
        end

        best_so_far(t) = bestFitness;
    end

    convergenceCurve = best_so_far;

end

%% --- Local Levy Function ---
function z = levy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta);
    u = random('Normal',0,sigma_u,n,m);
    v = random('Normal',0,1,n,m);
    z = u ./ (abs(v).^(1/beta));
end
