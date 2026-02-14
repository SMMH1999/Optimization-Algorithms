function [bestFitness, bestPosition, convergenceCurve] = ECBSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Enhanced Connected Banking System Optimizer (ECBSO)
    % =========================================================================
    % Author: (Based on original ECBSO source code)
    %
    % Description:
    % ECBSO is a hybrid population-based metaheuristic that combines
    % Connected Banking System Optimizer (CBSO), Equilibrium Optimizer (EO),
    % Gaussian Local Search (GLS), Lévy flight, and Brownian motion mechanisms
    % to balance exploration and exploitation.
    %
    % Inputs:
    %   lb        : Lower bound (scalar or 1×dim vector)
    %   ub        : Upper bound (scalar or 1×dim vector)
    %   dim       : Number of decision variables (dimension)
    %   nPop      : Population size
    %   maxItr    : Maximum number of iterations
    %   objFun    : Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      : Best objective value found
    %   bestPosition     : Best solution vector found (1×dim)
    %   convergenceCurve : Best fitness value at each iteration (1×maxItr)
    %
    % Tunable Parameters:
    %   A       : GLS scaling factor
    %   Cmax    : GLS memory size
    %   a1, a2  : EO control parameters
    %   GP      : EO generation probability
    % =========================================================================

    %% ============================ Initialization ============================

    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    LB = repmat(lb, nPop, 1);
    UB = repmat(ub, nPop, 1);

    X = rand(nPop, dim) .* (UB - LB) + LB;
    fitness = inf(nPop, 1);

    bestPosition = zeros(1, dim);
    bestFitness = inf;
    convergenceCurve = zeros(1, maxItr);

    % GLS parameters
    lu = [lb; ub];
    A = 30;
    C = 0;
    Cmax = 3000;
    St = zeros(Cmax, dim);
    B = 200 ./ (ub - lb);

    % EO parameters
    V = 1;  a1 = 2;  a2 = 1;  GP = 0.5;
    Ceq2 = zeros(1,dim);  Ceq2_fit = inf;
    Ceq3 = zeros(1,dim);  Ceq3_fit = inf;
    Ceq4 = zeros(1,dim);  Ceq4_fit = inf;

    %% ========================== Fitness Evaluation ==========================

    for i = 1:nPop
        fitness(i) = objFun(X(i,:));
        if fitness(i) < bestFitness
            bestFitness = fitness(i);
            bestPosition = X(i,:);
        end
    end

    fit_old = fitness;
    X_old = X;

    %% ============================= Main Loop ================================

    for t = 1:maxItr

        %% --------- Ranking and Equilibrium Pool Update ---------
        for i = 1:nPop
            if fitness(i) < bestFitness
                bestFitness = fitness(i);  bestPosition = X(i,:);
            elseif fitness(i) < Ceq2_fit
                Ceq2_fit = fitness(i);  Ceq2 = X(i,:);
            elseif fitness(i) < Ceq3_fit
                Ceq3_fit = fitness(i);  Ceq3 = X(i,:);
            elseif fitness(i) < Ceq4_fit
                Ceq4_fit = fitness(i);  Ceq4 = X(i,:);
            end
        end

        Ceq_ave = (bestPosition + Ceq2 + Ceq3 + Ceq4) / 4;
        C_pool = [bestPosition; Ceq2; Ceq3; Ceq4; Ceq_ave];

        tt = (1 - t/maxItr)^(a2 * t/maxItr);
        lambda = rand(nPop, dim);
        r1 = rand(nPop,1);
        r2 = rand(nPop,1);
        rn = randi(size(C_pool,1), nPop, 1);
        rr = rand(nPop, dim);

        [~, index] = sort(fitness, 'descend');
        rank = zeros(1,nPop);
        for k = 1:nPop
            rank(index(k)) = (nPop-k+1)/nPop;
        end

        Elite = repmat(bestPosition, nPop, 1);
        CF = rand * (1 - t/maxItr);

        Levy_step = levy(nPop, dim, 1.5);
        Brown_step = randn(nPop, dim);

        %% ---------------- Population Update ----------------
        for i = 1:nPop
            if rank(i) > 0.5
                Ceq = C_pool(rn(i),:);
                F = a1 * sign(rr(i,:) - 0.5) .* (exp(-lambda(i,:) .* tt) - 1);
                GCP = 0.5 * r1(i) * ones(1,dim) * (r2(i) >= GP);
                G = (GCP .* (Ceq - lambda(i,:) .* X(i,:))) .* F;
                X(i,:) = Ceq + (X(i,:) - Ceq).*F + (G./lambda(i,:) * V).*(1 - F);
            else
                SEL1 = randi(nPop);
                SEL2 = randi(nPop);
                if t < 0.4*maxItr
                    temp_pos = bestPosition + randn(1,dim);
                end
                for j = 1:dim
                    if t < 0.2*maxItr
                        X(i,j) = temp_pos(j);
                    elseif t <= 0.4*maxItr
                        if i > nPop/2
                            X(i,j) = Elite(i,j) + CF*(Brown_step(i,j)*(rand*Elite(i,j)-X(SEL1,j)));
                        else
                            X(i,j) = temp_pos(j);
                        end
                    else
                        X(i,j) = Elite(i,j) + CF*(Levy_step(i,j)*(Levy_step(i,j)*X(SEL1,j)-X(SEL2,j)));
                    end
                end
            end
        end

        %% ---------------- Boundary Control ----------------
        X = max(X, LB);
        X = min(X, UB);

        %% ---------------- Fitness Update ----------------
        for i = 1:nPop
            fitness(i) = objFun(X(i,:));
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = X(i,:);
            end

            C = C + 1;
            St(mod(C-1,Cmax)+1,:) = X(i,:);
            if C >= Cmax
                V0 = std(St,0,1).*B;
                C = 0;
                X = X + randn(nPop,dim).*V0;
                X = max(X, LB);
                X = min(X, UB);
            end
        end

        %% --------- Elitism Memory ---------
        Inx = (fit_old < fitness);
        X(Inx,:) = X_old(Inx,:);
        fitness(Inx) = fit_old(Inx);
        fit_old = fitness;
        X_old = X;

        %% --------- Random Perturbation ---------
        if rand < 0.2
            X = X + rand*(CF*LB + rand(nPop,dim).*(UB-LB));
        else
            X = X + rand*rand*(X(randperm(nPop),:) - X(randperm(nPop),:));
        end

        %% ---------------- Convergence Curve ----------------
        convergenceCurve(t) = bestFitness;

    end

end

%% ============================= Levy Flight ==============================

function z = levy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta);
    u = randn(n,m)*sigma_u;
    v = randn(n,m);
    z = u ./ (abs(v).^(1/beta));
end
