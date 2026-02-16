function [bestFitness, bestPosition, convergenceCurve] = LSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Light Spectrum Optimizer (LSO)
    % =========================================================================
    % Authors: Reda Mohamed & Mohamed Abdel-Basset
    % Reference:
    % Abdel-Basset, M., Mohamed, R.
    % Light Spectrum Optimizer: A Novel Physics-Inspired Metaheuristic
    % Optimization Algorithm, Mathematics 2022, 10(19), 3466.
    %
    % -------------------------------------------------------------------------
    % Description:
    % LSO is a physics-inspired metaheuristic based on light dispersion,
    % refraction, reflection, and scattering phenomena. The algorithm models
    % colorful dispersion and scattering stages to balance exploration and
    % exploitation in the search space.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    % lb        : Lower bound (1×dim vector or scalar)
    % ub        : Upper bound (1×dim vector or scalar)
    % dim       : Number of decision variables
    % nPop      : Number of search agents (Light Rays)
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % -------------------------------------------------------------------------
    % Outputs:
    % bestFitness      : Best fitness value found
    % bestPosition     : Best solution vector found
    % convergenceCurve : Best fitness at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    % Ps  : Probability of first/second scattering stages
    % Pe  : Exchange control parameter between scattering strategies
    % Ph  : Hybrid boundary handling probability
    % B   : Exploitation probability in first scattering stage
    % =========================================================================

    %% Parameter Settings
    red = 1.3318;
    violet = 1.3435;

    Ps = 0.05;
    Pe = 0.6;
    Ph = 0.4;
    B  = 0.05;

    %% Bound Handling
    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    %% Initialization
    LightRays = rand(nPop, dim) .* (ub - lb) + lb;
    fitness   = zeros(nPop,1);

    for i = 1:nPop
        fitness(i) = objFun(LightRays(i,:));
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = LightRays(idx,:);

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    t = 1;
    while t <= maxItr
        for i = 1:nPop

            nA = LightRays(randi(nPop),:);
            nB = LightRays(i,:);
            nC = bestPosition;
            xbar = mean(LightRays,1);

            norm_nA = nA / norm(nA);
            norm_nB = nB / norm(nB);
            norm_nC = nC / norm(nC);
            Incid_norm = xbar / norm(xbar);

            k = red + rand*(violet-red);
            p = rand;
            q = rand;

            L1 = (1/k)*(Incid_norm - norm_nA.*dot(norm_nA,Incid_norm)) ...
                - norm_nA.*sqrt(abs(1-(1/k^2)+(1/k^2)*dot(norm_nA,Incid_norm)^2));

            L2 = L1 - 2*norm_nB.*dot(L1,norm_nB);

            L3 = k*(L2 - norm_nC.*dot(norm_nC,L2)) ...
                + norm_nC.*sqrt(abs(1-k^2 + k^2*dot(norm_nC,L2)^2));

            a = rand*(1 - t/maxItr);
            ginv = gammaincinv(a,1);
            GI = a*(1/rand)*ginv;
            Epsln = a*randn(1,dim);

            if p <= q
                New = LightRays(i,:) + GI.*Epsln.*rand(1,dim).*(L1-L3) ...
                    .*(LightRays(randi(nPop),:)-LightRays(randi(nPop),:));
            else
                New = LightRays(i,:) + GI.*Epsln.*rand(1,dim).*(L2-L3) ...
                    .*(LightRays(randi(nPop),:)-LightRays(randi(nPop),:));
            end

            % Boundary control
            if rand < Ph
                New = max(New, lb);
                New = min(New, ub);
            else
                for j = 1:dim
                    if New(j) > ub(j) || New(j) < lb(j)
                        New(j) = lb(j) + rand*(ub(j)-lb(j));
                    end
                end
            end

            Fnew = objFun(New);

            if Fnew <= fitness(i)
                LightRays(i,:) = New;
                fitness(i) = Fnew;
            end

            if Fnew <= bestFitness
                bestFitness = Fnew;
                bestPosition = New;
            end

            % -------- Scattering Stage --------
            sortedFit = sort(fitness);
            F = abs((fitness(i)-bestFitness) / ...
                (bestFitness - sortedFit(end) + eps));

            if F < rand || rand < Ps
                if rand < Pe
                    New = LightRays(i,:) + rand*(LightRays(randi(nPop),:) ...
                        - LightRays(randi(nPop),:)) ...
                        + (rand < B).*rand(1,dim).*(bestPosition-LightRays(i,:));
                else
                    New = (2*cos(rand*pi)).*(bestPosition .* LightRays(i,:));
                end
            else
                U = rand(1,dim) > rand(1,dim);
                New = U.*(LightRays(randi(nPop),:) + abs(randn) ...
                    .*(LightRays(randi(nPop),:) - LightRays(randi(nPop),:))) ...
                    + (~U).*LightRays(i,:);
            end

            % Boundary control
            if rand < Ph
                New = max(New, lb);
                New = min(New, ub);
            else
                for j = 1:dim
                    if New(j) > ub(j) || New(j) < lb(j)
                        New(j) = lb(j) + rand*(ub(j)-lb(j));
                    end
                end
            end

            Fnew = objFun(New);

            if Fnew <= fitness(i)
                LightRays(i,:) = New;
                fitness(i) = Fnew;
            end

            if Fnew <= bestFitness
                bestFitness = Fnew;
                bestPosition = New;
            end

            convergenceCurve(t) = bestFitness;
            t = t + 1;

            if t > maxItr
                break;
            end
        end
    end

end
