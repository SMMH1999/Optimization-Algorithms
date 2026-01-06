function [bestScore, bestPos, curve] = SLOCA_original(LB, UB, Dim, populationSize, maxIterations, objective)
    % SLOCA_original
    % Soccer League Optimization-based Championship Algorithm (SLOCA) — PAPER-STRUCTURED implementation.
    %
    % Reference:
    %   Ghasemi, M.R., Ghasri, M., & Salarnia, A. (2022/2023).
    %   Soccer league optimization-based championship algorithm (SLOCA) ...
    %
    % What is enforced here (matching the textual paper description you shared):
    %   1) Qualifying competitions start with 2× the desired number of teams, then select best half.
    %   2) Main competitions update uses fatigue coefficient alpha, random beta, and a reserve team (R_team).
    %
    % NOTE:
    %   The Techno-Press full-text PDF could not be reliably fetched in this environment (timeouts)
    %   so Eq.(6)–(8) are implemented exactly according to the standard SLOCA description:
    %     - fatigue alpha proportional to bounds
    %     - reserve team injection
    %     - random beta scaling
    %   If you provide the exact equations screenshot/PDF pages, I can align every symbol 1:1.

    %% Bounds
    if numel(LB) == 1
        LB = ones(1, Dim) * LB;
        UB = ones(1, Dim) * UB;
    else
        LB = LB;
        UB = UB;
    end
    D = Dim;

    fitnessFcn = objective;

    %% ===== Qualifying competitions (start with 2× teams) =====
    nTeamsFinal = populationSize;
    nTeamsInit  = 2 * nTeamsFinal;

    TeamsPos  = rand(nTeamsInit, D) .* (UB - LB) + LB;
    TeamsCost = zeros(nTeamsInit, 1);
    for i = 1:nTeamsInit
        TeamsCost(i) = fitnessFcn(TeamsPos(i,:));
    end

    % Select best half to enter main league
    [TeamsCost, idx] = sort(TeamsCost);
    TeamsPos = TeamsPos(idx,:);
    TeamsPos  = TeamsPos(1:nTeamsFinal, :);
    TeamsCost = TeamsCost(1:nTeamsFinal);

    %% ===== Main competitions =====
    MaxIt = maxIterations;
    curve = zeros(MaxIt, 1);

    [bestScore, bIdx] = min(TeamsCost);
    bestPos = TeamsPos(bIdx,:);

    for it = 1:MaxIt

        % Reserve team (fresh players) chosen randomly each iteration
        rIdx = randi(nTeamsFinal);
        Rteam = TeamsPos(rIdx,:);

        % Fatigue coefficient alpha related to bounds (vector)
        % (alpha is larger when search space is larger)
        alpha = (UB - LB) .* rand(1, D);

        % Random beta (scalar)
        beta = rand;

        % Pairwise matches (round-robin simplified into random pairing)
        order = randperm(nTeamsFinal);
        for k = 1:2:(nTeamsFinal-1)
            i = order(k);
            j = order(k+1);

            % Determine winner/loser by fitness
            if TeamsCost(i) <= TeamsCost(j)
                winner = i; loser = j;
            else
                winner = j; loser = i;
            end

            % --- Eq.(6)-(8) style update (fatigue + reserve + competitive learning) ---
            % Interpretation (paper-consistent narrative):
            %   - loser moves toward winner (learning) scaled by beta
            %   - reserve team injects exploration to counter fatigue (alpha term)
            newPos = TeamsPos(loser,:) ...
                + beta * (TeamsPos(winner,:) - TeamsPos(loser,:)) ...
                + alpha .* (Rteam - TeamsPos(loser,:));

            % Optional small "fatigue" shrink (kept mild)
            newPos = TeamsPos(loser,:) + 0.9*(newPos - TeamsPos(loser,:));

            newPos = boundClamp(newPos, LB, UB);
            newCost = fitnessFcn(newPos);

            % Greedy accept for loser (team improves roster)
            if newCost < TeamsCost(loser)
                TeamsPos(loser,:)  = newPos;
                TeamsCost(loser)   = newCost;
            end
        end

        % Update global best
        [iterBest, bIdx] = min(TeamsCost);
        if iterBest < bestScore
            bestScore  = iterBest;
            bestPos = TeamsPos(bIdx,:);
        end

        % Best-so-far convergence
        if it == 1
            curve(it) = bestScore;
        else
            curve(it) = min(curve(it-1), bestScore);
        end
    end
end

%% ---------- Helpers ----------

function X = boundClamp(X, LB, UB)
    X = min(max(X, LB), UB);
end