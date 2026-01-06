function [bestScore, bestPos, curve] = WHO(LB, UB, Dim, populationNo, maxItr, objective)
% WHO - Wild Horse Optimizer (refactored, structure-fixed)
%
% Signature (framework-compatible):
%   [bestScore, bestPos, curve] = WHO(LB, UB, Dim, populationNo, maxItr, objective)
%
% Notes:
% - Core mechanics preserved: Stallions, groups (foals), TDR-based update, crossover, exchange.
% - Cleaned grouping + fixed r3/rr scope bug in stallion update.

    %% ---- Bounds normalize ----
    if isscalar(LB), LB = repmat(LB, 1, Dim); else, LB = LB(:).'; end
    if isscalar(UB), UB = repmat(UB, 1, Dim); else, UB = UB(:).'; end

    %% ---- Parameters (original-style constants) ----
    PS = 0.2;   % Stallions Percentage
    PC = 0.13;  % Crossover Probability

    nStallion = max(1, ceil(PS * populationNo));
    nFoal     = max(0, populationNo - nStallion);

    %% ---- Data templates ----
    AgentTemplate.pos  = zeros(1, Dim);
    AgentTemplate.cost = inf;

    StallionTemplate.pos  = zeros(1, Dim);
    StallionTemplate.cost = inf;
    StallionTemplate.group = repmat(AgentTemplate, 0, 1);

    %% ---- Initialize Foals ----
    foals = repmat(AgentTemplate, nFoal, 1);
    for i = 1:nFoal
        foals(i).pos  = randPos(LB, UB, Dim);
        foals(i).cost = objective(foals(i).pos);
    end

    %% ---- Initialize Stallions ----
    stallions = repmat(StallionTemplate, nStallion, 1);
    for i = 1:nStallion
        stallions(i).pos  = randPos(LB, UB, Dim);
        stallions(i).cost = objective(stallions(i).pos);
    end

    %% ---- Assign foals to stallions (round-robin after shuffle) ----
    if nFoal > 0
        foals = foals(randperm(nFoal));
        for j = 1:nFoal
            s = mod(j-1, nStallion) + 1;
            stallions(s).group(end+1, 1) = foals(j); %#ok<AGROW>
        end
    end

    %% ---- Exchange (ensure each stallion is best in its own group) ----
    stallions = exchange(stallions);

    %% ---- Global best (best stallion) ----
    [bestScore, bestIdx] = min([stallions.cost]);
    bestPos = stallions(bestIdx).pos;

    curve = inf(maxItr, 1);
    curve(1) = bestScore;

    %% ===================== Main loop =====================
    for it = 2:maxItr
        TDR = 1 - it / maxItr;  % decreases over time

        % cache current global best for this iteration
        globalLeader.pos  = bestPos;
        globalLeader.cost = bestScore;

        for s = 1:nStallion
            % ---- Sort group (best to worst) ----
            if ~isempty(stallions(s).group)
                [~, ord] = sort([stallions(s).group.cost], 'ascend');
                stallions(s).group = stallions(s).group(ord);
            end

            ng = numel(stallions(s).group);

            % ---- Update each foal in the group ----
            for j = 1:ng
                if rand > PC
                    % main update (TDR-controlled)
                    [r3, rr] = mixedRandVector(TDR, Dim);

                    newPos = 2 .* r3 .* cos(2*pi.*rr) .* ...
                             (stallions(s).pos - stallions(s).group(j).pos) + stallions(s).pos;
                else
                    % crossover: average of worst foals from two other stallions
                    if nStallion >= 3
                        idx = randperm(nStallion);
                        idx(idx == s) = [];
                        a = idx(1);
                        c = idx(2);

                        x1 = worstMemberPos(stallions(c), Dim, LB, UB);
                        x2 = worstMemberPos(stallions(a), Dim, LB, UB);
                        newPos = (x1 + x2) / 2;
                    else
                        % fallback when not enough stallions
                        newPos = randPos(LB, UB, Dim);
                    end
                end

                newPos = clamp(newPos, LB, UB);
                stallions(s).group(j).pos  = newPos;
                stallions(s).group(j).cost = objective(newPos);
            end

            % ---- Update stallion position toward/around global best ----
            % (Bug fix: r3/rr must be defined here, not “leaked” from foal loop)
            [r3s, rrs] = mixedRandVector(TDR, Dim);

            if rand < 0.5
                candidate = 2 .* r3s .* cos(2*pi.*rrs) .* (globalLeader.pos - stallions(s).pos) + globalLeader.pos;
            else
                candidate = 2 .* r3s .* cos(2*pi.*rrs) .* (globalLeader.pos - stallions(s).pos) - globalLeader.pos;
            end

            candidate = clamp(candidate, LB, UB);
            candCost  = objective(candidate);

            if candCost < stallions(s).cost
                stallions(s).pos  = candidate;
                stallions(s).cost = candCost;
            end
        end

        % ---- Exchange again (leader within each group) ----
        stallions = exchange(stallions);

        % ---- Update global best ----
        [iterBest, bestIdx] = min([stallions.cost]);
        if iterBest < bestScore
            bestScore = iterBest;
            bestPos   = stallions(bestIdx).pos;
        end

        curve(it) = min(curve(it-1), bestScore);
    end
end

%% ===================== Helpers =====================

function X = randPos(LB, UB, D)
    X = LB + rand(1, D) .* (UB - LB);
end

function X = clamp(X, LB, UB)
    X = min(max(X, LB), UB);
end

function [r3, rr] = mixedRandVector(TDR, D)
% Build r3 vector exactly in the original spirit:
%   z ~ Bernoulli(TDR), r1 scalar, r2 vector, choose componentwise.
% Then rr = -2 + 4*r3 (vector)

    z  = rand(1, D) < TDR;
    r1 = rand;         % scalar
    r2 = rand(1, D);   % vector

    r3 = r2;
    r3(~z) = r1;
    rr = -2 + 4 .* r3;
end

function x = worstMemberPos(stallion, D, LB, UB)
% If group exists, return the worst (last after sorting asc); else random.
    if isempty(stallion.group)
        x = LB + rand(1, D) .* (UB - LB);
    else
        % assume group already sorted in many places; safe to compute anyway
        [~, ord] = sort([stallion.group.cost], 'ascend');
        g = stallion.group(ord);
        x = g(end).pos;
    end
end

function stallions = exchange(stallions)
% Swap stallion with best foal in its group if the foal is better
    for i = 1:numel(stallions)
        if isempty(stallions(i).group), continue; end

        [bestFoalCost, idx] = min([stallions(i).group.cost]);
        if bestFoalCost < stallions(i).cost
            bestFoal = stallions(i).group(idx);

            % put stallion into that foal slot
            stallions(i).group(idx).pos  = stallions(i).pos;
            stallions(i).group(idx).cost = stallions(i).cost;

            % stallion becomes the best foal
            stallions(i).pos  = bestFoal.pos;
            stallions(i).cost = bestFoal.cost;
        end
    end
end
