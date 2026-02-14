function [bestFitness, bestPosition, convergenceCurve] = EVO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Energy Valley Optimizer (EVO)
    % =========================================================================
    % Algorithm: Energy Valley Optimizer (EVO)
    % Reference: Azizi M., Aickelin U., Khorshidi H., Baghalzadeh M.
    % "Energy valley optimizer: a novel metaheuristic algorithm for global
    % and engineering optimization." Scientific Reports, 13, 226 (2023).
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size (number of particles)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Description:
    % EVO is a population-based metaheuristic inspired by energy valley
    % mechanisms. Each candidate solution updates its position based on
    % relative energy levels, neighborhood information, and the global best.
    % The algorithm alternates between exploitation toward the best solution
    % and exploration via stochastic perturbations.
    %
    % -------------------------------------------------------------------------
    % =========================================================================

    %% =========================  Bounds Handling  ============================
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% =========================  Initialization  =============================
    Particles = rand(nPop, dim) .* (ub - lb) + lb;
    NELs = zeros(nPop,1);

    for i = 1:nPop
        NELs(i) = objFun(Particles(i,:));
    end

    % Sort population
    [NELs, sortIdx] = sort(NELs);
    Particles = Particles(sortIdx,:);

    BS = Particles(1,:);          % Best Solution
    BS_NEL = NELs(1);             % Best Fitness
    WS_NEL = NELs(end);           % Worst Fitness

    convergenceCurve = zeros(1, maxItr);

    %% =========================  Main Loop  ==================================
    for t = 1:maxItr

        NewParticles = [];
        NewNELs = [];

        for i = 1:nPop

            % --- Distance Calculation ---
            Dist = sqrt(sum((Particles - Particles(i,:)).^2, 2));
            [~, a] = sort(Dist);

            % --- Neighbor Selection ---
            CnPtIndex = randi(nPop);
            if CnPtIndex < 3
                CnPtIndex = CnPtIndex + 2;
            end

            CnPtA = Particles(a(2:CnPtIndex),:);
            X_NG = mean(CnPtA, 1);
            X_CP = mean(Particles, 1);
            EB = mean(NELs);

            SL = (NELs(i) - BS_NEL) / (WS_NEL - BS_NEL + eps);
            SB = rand;

            if NELs(i) > EB

                if SB > SL

                    % --- Exploitation toward BS and NG ---
                    AlphaIndex1 = randi(dim);
                    AlphaIndex2 = randi([1 dim], AlphaIndex1 , 1);
                    NP1 = Particles(i,:);
                    NP1(AlphaIndex2) = BS(AlphaIndex2);

                    GamaIndex1 = randi(dim);
                    GamaIndex2 = randi([1 dim], GamaIndex1 , 1);
                    NP2 = Particles(i,:);
                    NP2(GamaIndex2) = X_NG(GamaIndex2);

                    NP = [NP1; NP2];

                else

                    % --- Directed Movement ---
                    Ir = rand(1,2);
                    Jr = rand(1,dim);
                    NP1 = Particles(i,:) + (Jr .* (Ir(1)*BS - Ir(2)*X_CP) / (SL + eps));

                    Ir = rand(1,2);
                    Jr = rand(1,dim);
                    NP2 = Particles(i,:) + (Jr .* (Ir(1)*BS - Ir(2)*X_NG));

                    NP = [NP1; NP2];

                end

            else

                % --- Random Exploration ---
                NP1 = Particles(i,:) + randn * SL .* (rand(1,dim).*(ub-lb)+lb);
                NP = NP1;

            end

            % --- Boundary Control ---
            NP = max(NP, lb);
            NP = min(NP, ub);

            % --- Fitness Evaluation ---
            for k = 1:size(NP,1)
                NewParticles = [NewParticles; NP(k,:)];
                NewNELs = [NewNELs; objFun(NP(k,:))];
            end

        end

        % --- Elitist Selection ---
        NewParticles = [NewParticles; Particles];
        NewNELs = [NewNELs; NELs];

        [NewNELs, sortIdx] = sort(NewNELs);
        NewParticles = NewParticles(sortIdx,:);

        Particles = NewParticles(1:nPop,:);
        NELs = NewNELs(1:nPop);

        BS = Particles(1,:);
        BS_NEL = NELs(1);
        WS_NEL = NELs(end);

        bestFitness = BS_NEL;
        bestPosition = BS;
        convergenceCurve(t) = bestFitness;

    end

end
