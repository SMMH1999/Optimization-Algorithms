function [bestFitness, bestPosition, convergenceCurve] = EDMO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Enhanced Dwarf Mongoose Optimization Algorithm (EDMO)
    % Abbreviation: EDMO
    %
    % Author: (Based on DMOA with OLOBL + DQTOS + BDFSS enhancements)
    %
    % Description:
    %   EDMO is an enhanced version of the Dwarf Mongoose Optimization Algorithm
    %   integrating three improvement strategies:
    %       1) OLOBL  - Orthogonal Learning Opposition-Based Learning
    %       2) DQTOS  - Dynamic Quantum Tunneling Optimization Strategy
    %       3) BDFSS  - Bionic Dynamic Focusing Search Strategy
    %
    % Input Parameters:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Output Parameters:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness at each iteration (maxItr×1)
    %
    % Tunable Internal Parameters:
    %   nBabysitter  - Number of babysitters (default = 3)
    %   peep         - Alpha female vocalization coefficient (default = 2)
    %   L            - Abandonment threshold for babysitters
    %   V0           - Initial potential barrier (DQTOS)
    %   tau_q        - Temperature decay constant
    %   alpha_min/max- BDFSS light response bounds
    %
    % =========================================================================

    %% ========================= Initialization ================================
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    VarSize = [1 dim];

    nBabysitter = 3;
    nAlphaGroup = nPop - nBabysitter;
    nScout = nAlphaGroup;
    L = round(0.6 * dim * nBabysitter);
    peep = 2;

    pop.Position = [];
    pop.Cost = [];

    pop = repmat(pop, nAlphaGroup, 1);

    BestSol.Cost = inf;
    tau = inf;
    sm = inf(nAlphaGroup, 1);

    for i = 1:nAlphaGroup
        pop(i).Position = lb + rand(VarSize) .* (ub - lb);
        pop(i).Position = min(max(pop(i).Position, lb), ub);
        pop(i).Cost = objFun(pop(i).Position);
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end

    C = zeros(nAlphaGroup,1);
    convergenceCurve = zeros(maxItr,1);

    %% =========================== Main Loop ===================================
    for t = 1:maxItr

        %% ---------------- Alpha Group Phase ----------------------------------
        costs = [pop.Cost]';
        MeanCost = mean(costs);
        F = exp(-costs/(MeanCost + eps));
        P = F / sum(F);

        for m = 1:nAlphaGroup
            i = RouletteWheelSelection(P);
            K = [1:i-1 i+1:nAlphaGroup];
            k = K(randi(numel(K)));

            phi = (peep/2) * (-1 + 2*rand(VarSize));
            newPos = pop(i).Position + phi .* (pop(i).Position - pop(k).Position);
            newPos = min(max(newPos, lb), ub);
            newCost = objFun(newPos);

            if newCost <= pop(i).Cost
                pop(i).Position = newPos;
                pop(i).Cost = newCost;
            else
                C(i) = C(i) + 1;
            end
        end

        %% ---------------- Scout Phase ----------------------------------------
        for i = 1:nScout
            K = [1:i-1 i+1:nAlphaGroup];
            k = K(randi(numel(K)));

            phi = (peep/2) * (-1 + 2*rand(VarSize));
            newPos = pop(i).Position + phi .* (pop(i).Position - pop(k).Position);
            newPos = min(max(newPos, lb), ub);
            newCost = objFun(newPos);

            sm(i) = (newCost - pop(i).Cost) / max(newCost, pop(i).Cost);

            if newCost <= pop(i).Cost
                pop(i).Position = newPos;
                pop(i).Cost = newCost;
            else
                C(i) = C(i) + 1;
            end
        end

        %% ---------------- Babysitter Replacement -----------------------------
        for i = 1:nBabysitter
            if C(i) >= L
                pop(i).Position = lb + rand(VarSize) .* (ub - lb);
                pop(i).Position = min(max(pop(i).Position, lb), ub);
                pop(i).Cost = objFun(pop(i).Position);
                C(i) = 0;
            end
        end

        %% ---------------- Global Best Update ---------------------------------
        for i = 1:nAlphaGroup
            if pop(i).Cost <= BestSol.Cost
                BestSol = pop(i);
            end
        end

        %% ---------------- Direction Control (Original DMO) -------------------
        CF = (1 - t/maxItr)^(2*t/maxItr);
        validSM = sm(~isinf(sm));
        newtau = mean(validSM);
        if isempty(newtau) || isnan(newtau)
            newtau = 0;
        end

        for i = 1:nScout
            phi = (peep/2) * (-1 + 2*rand(VarSize));
            M = (pop(i).Position .* sm(i)) ./ (pop(i).Position + eps);

            if newtau > tau
                cand = pop(i).Position - CF .* phi .* rand(VarSize) .* (pop(i).Position - M);
            else
                cand = pop(i).Position + CF .* phi .* rand(VarSize) .* (pop(i).Position - M);
            end

            tau = newtau;
            cand = min(max(cand, lb), ub);
            candCost = objFun(cand);

            if candCost <= pop(i).Cost
                pop(i).Position = cand;
                pop(i).Cost = candCost;
            end
        end

        %% ======================= OLOBL (Leader Only) =========================
        [olPos, olCost] = OLOBL(BestSol.Position, lb, ub, dim, t, maxItr, objFun);
        if olCost < BestSol.Cost
            BestSol.Position = olPos;
            BestSol.Cost = olCost;
            [~, worstIdx] = max([pop.Cost]);
            pop(worstIdx) = BestSol;
        end

        %% ======================= BDFSS Phase =================================
        costs = [pop.Cost];
        fbest = min(costs);
        fworst = max(costs);
        if fbest == fworst, fworst = fbest + 1e-12; end

        I = (fworst - costs) ./ (fworst - fbest + eps);
        alpha_min = 0.1; alpha_max = 0.9;
        alphas = alpha_min + (alpha_max-alpha_min) .* I;

        delta = 0.01 * (ub - lb);

        for i = 1:nScout
            grad = zeros(VarSize);
            for d = 1:dim
                e = zeros(VarSize); e(d) = 1;
                xp = min(max(pop(i).Position + delta.*e, lb), ub);
                xm = min(max(pop(i).Position - delta.*e, lb), ub);
                Ip = (fworst - objFun(xp)) / (fworst - fbest + eps);
                Im = (fworst - objFun(xm)) / (fworst - fbest + eps);
                grad(d) = (Ip - Im) / (2*delta(d));
            end

            social = (BestSol.Position - pop(i).Position);
            xi = rand();
            phi_local = rand(VarSize);
            cand = pop(i).Position + phi_local .* (social + xi .* alphas(i) .* grad);
            cand = min(max(cand, lb), ub);
            candCost = objFun(cand);

            if candCost < pop(i).Cost
                pop(i).Position = cand;
                pop(i).Cost = candCost;
                if candCost < BestSol.Cost
                    BestSol = pop(i);
                end
            end
        end

        %% ======================= DQTOS Phase =================================
        costs = [pop.Cost];
        fbest = min(costs);
        fworst = max(costs);

        V0 = 0.5;
        tau_q = 0.1 * maxItr;
        Tt = 1 / (1 + exp(t/tau_q));
        Delta = (ub - lb)/2 .* (1 - t/maxItr)^2;
        beta = 1;

        for i = 1:nAlphaGroup
            Vi = abs(pop(i).Cost - fbest) / (fworst - fbest + eps) * V0;
            Ek = Tt + eps;
            Pt = exp(-max(0, Vi/Ek));

            if rand < Pt
                dir = 0.7*randn(VarSize) + 0.3*(BestSol.Position - pop(i).Position);
                cand = pop(i).Position + beta*Tt.*Delta.*dir;
                cand = min(max(cand, lb), ub);
                candCost = objFun(cand);

                if candCost < pop(i).Cost
                    pop(i).Position = cand;
                    pop(i).Cost = candCost;
                    if candCost < BestSol.Cost
                        BestSol = pop(i);
                    end
                end
            end
        end

        %% ---------------- Convergence Recording -------------------------------
        convergenceCurve(t) = BestSol.Cost;
    end

    bestFitness = BestSol.Cost;
    bestPosition = BestSol.Position;

end

%% ======================= Local Functions =================================
function i = RouletteWheelSelection(P)
    r = rand;
    C = cumsum(P);
    i = find(r <= C, 1, 'first');
end

function [bestPos,bestCost] = OLOBL(x,lb,ub,D,t,Tmax,objFun)
    M = max(4, 2*ceil(log2(D+1)));
    mid = (ub + lb)/2;
    k = (1 + (t/max(1,Tmax))^2)^10;
    x_reflect = mid + mid./k - x./k;
    x_reflect = min(max(x_reflect, lb), ub);

    candCost = inf(M,1);
    candPos = zeros(M,D);

    for m = 1:M
        mask = rand(1,D) > 0.5;
        xp = x;
        xp(mask) = x_reflect(mask);
        xp = min(max(xp, lb), ub);
        candPos(m,:) = xp;
        candCost(m) = objFun(xp);
    end

    origCost = objFun(x);
    [bestCost, idx] = min([origCost; candCost]);
    if idx == 1
        bestPos = x;
    else
        bestPos = candPos(idx-1,:);
    end
end
