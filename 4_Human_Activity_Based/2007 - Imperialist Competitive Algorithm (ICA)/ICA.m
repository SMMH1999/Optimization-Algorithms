function [bestScore, bestPos, curve] = ICA(LB, UB, Dim, populationSize, maxIterations, objective)
    % ICA_original
    % Imperialist Competitive Algorithm (ICA) — ORIGINAL (paper-faithful) implementation.
    %
    % Reference:
    %   Atashpaz-Gargari, E., & Lucas, C. (2007). Imperialist Competitive Algorithm:
    %   An algorithm for optimization inspired by imperialistic competition.
    %
    % Key operators in original ICA include Assimilation, Revolution, and Imperialistic Competition. citeturn0search1
    %
    % IMPORTANT NOTES (for paper fidelity vs your previous refactor):
    %   - Assimilation uses a SCALAR distance x ~ U(0, beta*d) and angular deviation theta.
    %     Here we implement deviation by rotating the movement direction in a random plane.
    %   - Zeta (ξ) is set to 0.1 as commonly recommended in ICA references.
    %   - Revolution is included as a core operator (as described in ICA literature). citeturn0search1
    %
    % I/O signature kept compatible with your framework.
    % evaluateFitness(...) is only to support your mixed benchmark-call styles.

    %% Bounds
    if numel(LB) == 1
        VarMin = repmat(LB, 1, Dim);
        VarMax = repmat(UB, 1, Dim);
    else
        VarMin = LB;
        VarMax = UB;
    end
    D = Dim;

    %% Parameters (typical defaults used in original ICA implementations)
    Ncountry   = populationSize;
    Nimp       = max(round(0.1*Ncountry), 2);
    Ncol       = Ncountry - Nimp;
    MaxDecades = maxIterations;

    beta  = 2;        % assimilation coefficient (β)
    gamma = pi/4;     % deviation range: theta ~ U(-gamma, gamma)
    zeta  = 0.1;      % ξ (a.k.a zeta) for total cost

    RevolutionRate = 0.1;  % fraction of colonies revolved per decade (typical)
    UniteThresh    = 0.02; % fraction of search-space diagonal

    fitnessFcn = objective;

    %% Initialize countries
    Countries = rand(Ncountry, D) .* (VarMax - VarMin) + VarMin;
    Costs     = zeros(Ncountry, 1);
    for i = 1:Ncountry
        Costs(i) = fitnessFcn(Countries(i,:));
    end

    % Sort countries
    [Costs, idx] = sort(Costs);
    Countries    = Countries(idx,:);

    ImperialistsPos  = Countries(1:Nimp,:);
    ImperialistsCost = Costs(1:Nimp);

    ColoniesPos  = Countries(Nimp+1:end,:);
    ColoniesCost = Costs(Nimp+1:end);

    %% Create empires (allocate colonies based on normalized power)
    if max(ImperialistsCost) > 0
        Power = 1.3*max(ImperialistsCost) - ImperialistsCost;
    else
        Power = 0.7*max(ImperialistsCost) - ImperialistsCost;
    end
    Power = Power / sum(Power);

    NcolPerImp = round(Power * Ncol);
    % fix rounding
    NcolPerImp(end) = Ncol - sum(NcolPerImp(1:end-1));

    Empires = repmat(struct( ...
        'ImperialistPosition', [], 'ImperialistCost', [], ...
        'ColoniesPosition', [], 'ColoniesCost', [], ...
        'TotalCost', []), Nimp, 1);

    perm = randperm(Ncol);
    startIdx = 1;

    for i = 1:Nimp
        Empires(i).ImperialistPosition = ImperialistsPos(i,:);
        Empires(i).ImperialistCost     = ImperialistsCost(i);

        nci = NcolPerImp(i);
        if nci > 0
            colsIdx = perm(startIdx:startIdx+nci-1);
            startIdx = startIdx + nci;

            Empires(i).ColoniesPosition = ColoniesPos(colsIdx,:);
            Empires(i).ColoniesCost     = ColoniesCost(colsIdx,:);
        else
            Empires(i).ColoniesPosition = zeros(0, D);
            Empires(i).ColoniesCost     = zeros(0, 1);
        end

        Empires(i).TotalCost = Empires(i).ImperialistCost + zeta * mean_or_inf(Empires(i).ColoniesCost);
    end

    %% Initialize best-so-far
    [bestScore, bIdx] = min([Empires.ImperialistCost]);
    bestPos = Empires(bIdx).ImperialistPosition;

    curve = zeros(MaxDecades, 1);

    %% Main loop
    diagSpace = norm(VarMax - VarMin);
    Dth = UniteThresh * diagSpace;

    for dec = 1:MaxDecades

        %% Assimilation + Revolution + Possession + Total cost
        for i = 1:numel(Empires)
            % Assimilate colonies toward imperialist with angular deviation
            Empires(i) = AssimilateWithDeviation(Empires(i), beta, gamma, VarMin, VarMax);

            % Revolution: random reinitialization of some colonies
            Empires(i) = Revolution(Empires(i), RevolutionRate, VarMin, VarMax);

            % Re-evaluate colonies
            nC = size(Empires(i).ColoniesPosition, 1);
            Empires(i).ColoniesCost = zeros(nC, 1);
            for j = 1:nC
                Empires(i).ColoniesCost(j) = fitnessFcn(Empires(i).ColoniesPosition(j,:));
            end

            % Possession: swap best colony with imperialist if better
            Empires(i) = Possession(Empires(i));

            % Update total cost
            Empires(i).TotalCost = Empires(i).ImperialistCost + zeta * mean_or_inf(Empires(i).ColoniesCost);
        end

        %% Unite similar empires (distance between imperialists)
        Empires = UniteEmpires(Empires, Dth, zeta);

        %% Imperialistic competition: weakest empire loses a colony
        Empires = ImperialisticCompetition(Empires);

        %% Update global best
        [currentBest, bIdx] = min([Empires.ImperialistCost]);
        if currentBest < bestScore
            bestScore  = currentBest;
            bestPos = Empires(bIdx).ImperialistPosition;
        end

        %% Monotonic convergence
        if dec == 1
            curve(dec) = bestScore;
        else
            curve(dec) = min(curve(dec-1), bestScore);
        end
    end
end

%% ---------------- Helpers ----------------

function Empire = AssimilateWithDeviation(Empire, beta, gamma, VarMin, VarMax)
    if isempty(Empire.ColoniesPosition)
        return;
    end

    Imp = Empire.ImperialistPosition;
    Col = Empire.ColoniesPosition;

    nCol = size(Col, 1);
    for k = 1:nCol
        vec = Imp - Col(k,:);
        d   = norm(vec);

        if d < eps
            continue;
        end

        % x ~ U(0, beta*d)
        x = rand * beta * d;

        % theta ~ U(-gamma, gamma)
        theta = (2*rand - 1) * gamma;

        u = vec / d; % unit direction toward imperialist

        % Build a random orthonormal vector v perpendicular to u
        v = randn(1, numel(u));
        v = v - dot(v, u) * u;
        nv = norm(v);
        if nv < eps
            v = null_vec(u);
            nv = norm(v);
        end
        v = v / nv;

        % Rotate in the (u,v) plane by theta
        dir = cos(theta)*u + sin(theta)*v;

        Col(k,:) = Col(k,:) + x * dir;
    end

    % Bound handling (clamp)
    Col = min(max(Col, VarMin), VarMax);
    Empire.ColoniesPosition = Col;
end

function Empire = Revolution(Empire, rate, VarMin, VarMax)
    nCol = size(Empire.ColoniesPosition, 1);
    if nCol == 0
        return;
    end
    nRev = round(rate * nCol);
    if nRev < 1
        return;
    end
    idx = randperm(nCol, nRev);
    D = numel(VarMin);
    Empire.ColoniesPosition(idx,:) = rand(nRev, D) .* (VarMax - VarMin) + VarMin;
end

function Empire = Possession(Empire)
    if isempty(Empire.ColoniesCost)
        return;
    end
    [minC, idx] = min(Empire.ColoniesCost);
    if minC < Empire.ImperialistCost
        % swap
        tmpPos  = Empire.ImperialistPosition;
        tmpCost = Empire.ImperialistCost;

        Empire.ImperialistPosition = Empire.ColoniesPosition(idx,:);
        Empire.ImperialistCost     = Empire.ColoniesCost(idx);

        Empire.ColoniesPosition(idx,:) = tmpPos;
        Empire.ColoniesCost(idx)       = tmpCost;
    end
end

function Empires = UniteEmpires(Empires, Dth, zeta)
    i = 1;
    while i <= numel(Empires)
        j = i + 1;
        while j <= numel(Empires)
            if norm(Empires(i).ImperialistPosition - Empires(j).ImperialistPosition) <= Dth
                % merge weaker into stronger
                if Empires(i).ImperialistCost <= Empires(j).ImperialistCost
                    best = i; worst = j;
                else
                    best = j; worst = i;
                end

                Empires(best).ColoniesPosition = [Empires(best).ColoniesPosition;
                    Empires(worst).ImperialistPosition;
                    Empires(worst).ColoniesPosition];
                Empires(best).ColoniesCost = [Empires(best).ColoniesCost;
                    Empires(worst).ImperialistCost;
                    Empires(worst).ColoniesCost];

                Empires(best).TotalCost = Empires(best).ImperialistCost + zeta * mean_or_inf(Empires(best).ColoniesCost);

                Empires(worst) = [];
                % restart merge scan (simple)
                i = 1;
                j = 2;
                continue;
            end
            j = j + 1;
        end
        i = i + 1;
    end
end

function Empires = ImperialisticCompetition(Empires)
    if numel(Empires) <= 1
        return;
    end

    TotalCost = [Empires.TotalCost];
    [~, weakest] = max(TotalCost);

    % Possession probability based on normalized power
    P = max(TotalCost) - TotalCost;
    if all(P == 0)
        P = ones(size(P));
    end
    P = P / sum(P);

    winner = RouletteWheel(P);

    % Transfer one colony from weakest to winner (or imperialist if none)
    if ~isempty(Empires(weakest).ColoniesPosition)
        k = randi(size(Empires(weakest).ColoniesPosition, 1));

        Empires(winner).ColoniesPosition = [Empires(winner).ColoniesPosition;
            Empires(weakest).ColoniesPosition(k,:)];
        Empires(winner).ColoniesCost = [Empires(winner).ColoniesCost;
            Empires(weakest).ColoniesCost(k)];

        Empires(weakest).ColoniesPosition(k,:) = [];
        Empires(weakest).ColoniesCost(k)       = [];
    else
        % weakest has no colonies: transfer imperialist itself (empire collapses)
        Empires(winner).ColoniesPosition = [Empires(winner).ColoniesPosition;
            Empires(weakest).ImperialistPosition];
        Empires(winner).ColoniesCost = [Empires(winner).ColoniesCost;
            Empires(weakest).ImperialistCost];

        Empires(weakest) = [];
        return;
    end

    % eliminate powerless empires (no colonies and relatively weak)
    if isempty(Empires(weakest).ColoniesPosition)
        Empires(winner).ColoniesPosition = [Empires(winner).ColoniesPosition;
            Empires(weakest).ImperialistPosition];
        Empires(winner).ColoniesCost = [Empires(winner).ColoniesCost;
            Empires(weakest).ImperialistCost];
        Empires(weakest) = [];
    end
end

function idx = RouletteWheel(P)
    r = rand;
    C = cumsum(P);
    idx = find(r <= C, 1, 'first');
    if isempty(idx)
        idx = numel(P);
    end
end

function m = mean_or_inf(v)
    if isempty(v)
        m = inf;
    else
        m = mean(v);
    end
end

function v = null_vec(u)
    % fallback to find a vector not parallel to u
    D = numel(u);
    e = zeros(1, D);
    [~, k] = max(abs(u));
    k2 = mod(k, D) + 1;
    e(k2) = 1;
    v = e - dot(e,u)*u;
end

