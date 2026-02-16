function [bestScore, bestPos, curve] = ED(LB, UB, Dim, populationSize, maxIterations, objective)
    % ED_original
    % Enterprise Development-inspired metaheuristic (ED) — ORIGINAL (paper-faithful) implementation.
    %
    % Reference (as cited in your refactored code):
    %   Truong, D.-N., & Chou, J.-S. "Metaheuristic algorithm inspired by enterprise development ..."
    %   Engineering Structures (2024) (Article: 118679).
    %
    % IMPORTANT NOTES (for strict paper fidelity):
    %   - Structure phase (Eq. 3): uses ONE scalar rand(-1,1) multiplying the full direction vector.
    %   - Technology phase (Eq. 5): uses scalar random factors per vector term (not per-Dim).
    %   - People phase (Eq. 6): changes ONE Dim using scalar rand(-1,1) (already scalar by nature).
    %
    % I/O signature is kept compatible with your framework.
    % The only "compatibility" helper included is evaluateFitness(...) to support both CEC/ProbInfo
    % and (x') styles. It does NOT change the ED update equations.

    %% Bounds (scalar or vector)
    if numel(LB) == 1
        LB = ones(1, Dim) * LB;
        UB = ones(1, Dim) * UB;
    else
        LB = LB;
        UB = UB;
    end

    D   = Dim;
    NP  = populationSize;
    GEN = maxIterations;

    fitnessFcn = objective;

    %% Initialize population (Eq. 1)
    Pop    = zeros(NP, D);
    ObjVal = zeros(NP, 1);

    for i = 1:NP
        Pop(i, :) = LB + rand(1, D) .* (UB - LB);
        ObjVal(i) = fitnessFcn(Pop(i, :));
    end

    %% Initial best
    [GlobalMin, iBest] = min(ObjVal);
    Xbest = Pop(iBest, :);

    curve = zeros(GEN, 1);
    bestScore  = GlobalMin;
    bestPos = Xbest;

    %% Main loop
    for g = 1:GEN

        % Time varying parameter c(t) (as used in many ED pseudocodes)
        ct = 1 - rand * g / GEN;   %#ok<NASGU>

        if rand < 0.1
            % Task activity (Eq. 2): replace worst
            [~, idx] = sort(ObjVal);
            kWorst = idx(end);

            SolD = LB + rand(1, D) .* (UB - LB);
            f_d  = fitnessFcn(SolD);

            Pop(kWorst, :) = SolD;
            ObjVal(kWorst) = f_d;

            if f_d < GlobalMin
                GlobalMin = f_d;
                Xbest     = SolD;
            end

        else
            % Strategy selector a = ceil(3*c(t))
            a = ceil(3 * (1 - rand * g / GEN));  % keep the same "3-region" selection idea
            a = min(max(a, 1), 3);

            switch a
                case 1
                    %% Structure activity (Eq. 3) — SCALAR rand(-1,1)
                    for i = 1:NP
                        P = randperm(NP, 3);
                        h = P(1); p = P(2); k = P(3);

                        SolC_mean = (Pop(h, :) + Pop(p, :) + Pop(k, :)) / 3;

                        r = 2*rand - 1; % rand(-1,1) SCALAR
                        SolC = Pop(i, :) + r * (Xbest - SolC_mean);

                        SolC = boundRepair_resample(SolC, LB, UB);
                        f_c  = fitnessFcn(SolC);

                        if f_c <= GlobalMin
                            GlobalMin = f_c;
                            Xbest     = SolC;
                        end
                    end

                case 2
                    %% Technology activity (Eq. 5) — SCALAR random multipliers per term
                    for i = 1:NP
                        h = randi(NP);

                        r1 = rand;  % scalar
                        r2 = rand;  % scalar
                        SolB = Pop(i, :) + (r1 * (Xbest - Pop(i, :)) + r2 * (Xbest - Pop(h, :)));

                        SolB = boundRepair_resample(SolB, LB, UB);
                        f_b  = fitnessFcn(SolB);

                        if f_b <= ObjVal(i)
                            Pop(i, :) = SolB;
                            ObjVal(i) = f_b;

                            if f_b <= GlobalMin
                                GlobalMin = f_b;
                                Xbest     = SolB;
                            end
                        end
                    end

                case 3
                    %% People activity (Eq. 6)
                    for i = 1:NP
                        Change = randi(D);

                        A   = randperm(NP, 3);
                        nb1 = A(1); nb2 = A(2); nb3 = A(3);

                        SolA = Pop(i, :);

                        r = 2*rand - 1; % scalar rand(-1,1)
                        SolA(Change) = Pop(i, Change) + ...
                            (Xbest(Change) - (Pop(nb1, Change) + Pop(nb2, Change) + Pop(nb3, Change))/3) * r;

                        SolA = boundRepair_resample(SolA, LB, UB);
                        f_a  = fitnessFcn(SolA);

                        if f_a <= ObjVal(i)
                            Pop(i, :) = SolA;
                            ObjVal(i) = f_a;

                            if f_a <= GlobalMin
                                GlobalMin = f_a;
                                Xbest     = SolA;
                            end
                        end
                    end
            end
        end

        %% Track best-so-far (monotonic curve)
        if g == 1
            curve(g) = GlobalMin;
        else
            curve(g) = min(curve(g-1), GlobalMin);
        end

        bestScore  = curve(g);
        bestPos = Xbest;
    end
end

%% ---------- Helpers ----------

function X = boundRepair_resample(X, LB, UB)
    % Resample violated components uniformly within bounds (common in ED codes)
    lo = X < LB;
    hi = X > UB;

    if any(lo)
        X(lo) = LB(lo) + rand(1, nnz(lo)) .* (UB(lo) - LB(lo));
    end
    if any(hi)
        X(hi) = LB(hi) + rand(1, nnz(hi)) .* (UB(hi) - LB(hi));
    end
end
