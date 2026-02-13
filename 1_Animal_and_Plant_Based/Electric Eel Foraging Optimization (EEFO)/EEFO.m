function [bestFitness, bestPosition, convergenceCurve] = EEFO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Electric Eel Foraging Optimization (EEFO)
    % -------------------------------------------------------------------------
    % Author(s): W. Zhao, L. Wang, Z. Zhang, H. Fan, J. Zhang, S. Mirjalili,
    %            N. Khodadadi, Q. Cao
    % Source:
    % Electric eel foraging optimization: A new bio-inspired optimizer for
    % engineering applications, Expert Systems With Applications, 238 (2024),
    % 122200.
    %
    % Refactored to Benchmark Format (Single-Function Version)
    % -------------------------------------------------------------------------
    % INPUTS:
    %   lb        : Lower bound (scalar or 1×dim vector)
    %   ub        : Upper bound (scalar or 1×dim vector)
    %   dim       : Problem dimensionality (integer)
    %   nPop      : Population size (number of electric eels)
    %   maxItr    : Maximum number of iterations
    %   objFun    : Objective function handle (minimization)
    %
    % OUTPUTS:
    %   bestFitness      : Best objective value found
    %   bestPosition     : Best solution vector corresponding to bestFitness
    %   convergenceCurve : Best fitness value at each iteration (maxItr×1)
    %
    % -------------------------------------------------------------------------
    % Algorithm Parameters:
    %   E0     : Initial energy factor (time-varying)
    %   Alpha  : Resting area scale factor
    %   Beta   : Hunting area scale factor
    %   Eta    : Curling factor
    %
    % The algorithm simulates electric eel foraging behaviors including:
    %   - Interaction
    %   - Resting
    %   - Migrating
    %   - Hunting
    %
    % Minimization is assumed.
    % =========================================================================

    %% ========================= Initialization ==============================

    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';
        ub = ub(:)';
    end

    PopPos = rand(nPop, dim) .* (ub - lb) + lb;
    PopFit = zeros(1, nPop);

    for i = 1:nPop
        PopFit(i) = objFun(PopPos(i, :));
    end

    [bestFitness, idx] = min(PopFit);
    bestPosition = PopPos(idx, :);

    convergenceCurve = zeros(maxItr, 1);

    %% =========================== Main Loop =================================

    for t = 1:maxItr

        DirectVector = zeros(nPop, dim);
        E0 = 4 * sin(1 - t / maxItr);

        for i = 1:nPop

            E = E0 * log(1 / rand);  % Energy factor (Eq.30)

            % Direction Vector (Eq.6)
            if dim == 1
                DirectVector(i, :) = 1;
            else
                RandNum = ceil((maxItr - t) / maxItr * rand * (dim - 2) + 2);
                RandDim = randperm(dim);
                DirectVector(i, RandDim(1:RandNum)) = 1;
            end

            if E > 1
                % ===================== Interaction ===========================
                idxSet = [1:i-1, i+1:nPop];
                j = idxSet(randi(nPop - 1));

                if PopFit(j) < PopFit(i)
                    if rand > 0.5
                        newPos = PopPos(j,:) + randn .* DirectVector(i,:) .* ...
                            (mean(PopPos) - PopPos(i,:));
                    else
                        xr = rand(1, dim) .* (ub - lb) + lb;
                        newPos = PopPos(j,:) + randn .* DirectVector(i,:) .* ...
                            (xr - PopPos(i,:));
                    end
                else
                    if rand > 0.5
                        newPos = PopPos(i,:) + randn .* DirectVector(i,:) .* ...
                            (mean(PopPos) - PopPos(j,:));
                    else
                        xr = rand(1, dim) .* (ub - lb) + lb;
                        newPos = PopPos(i,:) + randn .* DirectVector(i,:) .* ...
                            (xr - PopPos(j,:));
                    end
                end

            else
                r = rand;

                if r < 1/3
                    % ======================= Resting ==========================
                    Alpha = 2 * (exp(1) - exp(t / maxItr)) * sin(2 * pi * rand);
                    rn = randi(nPop);
                    rd = randi(dim);

                    z = (PopPos(rn, rd) - lb(rd)) / (ub(rd) - lb(rd));
                    Z = lb + z .* (ub - lb);

                    Ri = Z + Alpha .* abs(Z - bestPosition);
                    newPos = Ri + randn * (Ri - round(rand) * PopPos(i,:));

                elseif r > 2/3
                    % ======================= Migrating ========================
                    rn = randi(nPop);
                    rd = randi(dim);

                    z = (PopPos(rn, rd) - lb(rd)) / (ub(rd) - lb(rd));
                    Z = lb + z .* (ub - lb);

                    Alpha = 2 * (exp(1) - exp(t / maxItr)) * sin(2 * pi * rand);
                    Ri = Z + Alpha .* abs(Z - bestPosition);

                    Beta = 2 * (exp(1) - exp(t / maxItr)) * sin(2 * pi * rand);
                    Hr = bestPosition + Beta .* abs(mean(PopPos) - bestPosition);

                    L = 0.01 * abs(levy(dim));
                    newPos = -rand * Ri + rand * Hr - L .* (Hr - PopPos(i,:));

                else
                    % ======================== Hunting =========================
                    Beta = 2 * (exp(1) - exp(t / maxItr)) * sin(2 * pi * rand);
                    Hprey = bestPosition + Beta .* abs(mean(PopPos) - bestPosition);

                    r4 = rand;
                    Eta = exp(r4 * (1 - t) / maxItr) * cos(2 * pi * r4);

                    newPos = Hprey + Eta * (Hprey - round(rand) * PopPos(i,:));
                end
            end

            % ==================== Boundary Control ==========================
            newPos = spaceBound(newPos, ub, lb);
            newFit = objFun(newPos);

            % ===================== Greedy Selection =========================
            if newFit < PopFit(i)
                PopFit(i) = newFit;
                PopPos(i,:) = newPos;

                if newFit <= bestFitness
                    bestFitness = newFit;
                    bestPosition = newPos;
                end
            end

        end

        % ================== Convergence Tracking ============================
        convergenceCurve(t) = bestFitness;

    end

end

%% ========================== Levy Flight ================================

function step = levy(d)
    b = 1.5;
    sigma = (gamma(1+b) * sin(pi*b/2) / ...
        (gamma((1+b)/2) * b * 2^((b-1)/2)))^(1/b);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    step = u ./ abs(v).^(1/b);
end

%% ========================= Boundary Handler ============================

function X = spaceBound(X, ub, lb)
    Dim = length(X);
    S = (X > ub) + (X < lb);
    X = (rand(1,Dim).*(ub-lb)+lb).*S + X.*(~S);
end
