function [bestFitness, bestPosition, convergenceCurve] = ASO(lb, ub, dim, nPop, maxItr, objFun)
    %--------------------------------------------------------------------------
    % Atom Search Optimization (ASO)
    %--------------------------------------------------------------------------
    % Developed based on:
    % W. Zhao, L. Wang and Z. Zhang, "Atom Search Optimization and its
    % Application to Solve a Hydrogeologic Parameter Estimation Problem",
    % Knowledge-Based Systems, 2019, 163:283-304.
    %--------------------------------------------------------------------------
    % Algorithm: Atom Search Optimization (ASO)
    %
    % Description:
    % ASO is a physics-inspired metaheuristic algorithm based on molecular
    % dynamics theory. Each candidate solution is treated as an atom whose
    % motion is governed by interaction forces (Lennard-Jones potential) and
    % attraction toward the best atom. The search process is guided by mass,
    % interaction force, and acceleration updates.
    %
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Dimension of decision variables
    %   nPop      - Number of atoms (population size)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (maxItr×1)
    %
    % Tunable Parameters (original ASO settings):
    %   alpha = 50       - Depth weight (interaction force factor)
    %   beta  = 0.2      - Multiplier weight (best attraction factor)
    %
    % Notes:
    %   - Minimization problem
    %   - Boundary handling via random reinitialization within bounds
    %   - Fully benchmark-compatible
    %--------------------------------------------------------------------------

    %% Parameters (original default values)
    alpha = 50;
    beta  = 0.2;

    %% Bound handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% Initialization
    Atom_Pop = rand(nPop, dim) .* (ub - lb) + lb;
    Atom_V   = rand(nPop, dim) .* (ub - lb) + lb;

    Fitness = zeros(nPop,1);
    for i = 1:nPop
        Fitness(i) = objFun(Atom_Pop(i,:));
    end

    convergenceCurve = zeros(maxItr,1);
    [bestFitness, idx] = min(Fitness);
    bestPosition = Atom_Pop(idx,:);
    convergenceCurve(1) = bestFitness;

    %% Initial Acceleration
    Atom_Acc = Acceleration(Atom_Pop, Fitness, 1, maxItr, dim, nPop, bestPosition, alpha, beta);

    %% Main Loop
    for t = 2:maxItr

        convergenceCurve(t) = convergenceCurve(t-1);

        % Velocity & Position Update
        Atom_V   = rand(nPop, dim) .* Atom_V + Atom_Acc;
        Atom_Pop = Atom_Pop + Atom_V;

        % Boundary Control + Fitness Evaluation
        for i = 1:nPop
            TU = Atom_Pop(i,:) > ub;
            TL = Atom_Pop(i,:) < lb;
            Atom_Pop(i,:) = Atom_Pop(i,:).*(~(TU+TL)) + ...
                (rand(1,dim).*(ub-lb)+lb).*(TU+TL);
            Fitness(i) = objFun(Atom_Pop(i,:));
        end

        [currentBest, idx] = min(Fitness);

        if currentBest < convergenceCurve(t)
            convergenceCurve(t) = currentBest;
            bestPosition = Atom_Pop(idx,:);
        else
            r = randi(nPop);
            Atom_Pop(r,:) = bestPosition;
        end

        % Update Acceleration
        Atom_Acc = Acceleration(Atom_Pop, Fitness, t, maxItr, dim, nPop, bestPosition, alpha, beta);
    end

    bestFitness = convergenceCurve(maxItr);

end

%% ------------------------------------------------------------------------
function Acc = Acceleration(Atom_Pop, Fitness, Iteration, Max_Iteration, Dim, Atom_Num, X_Best, alpha, beta)

    % Mass Calculation
    if max(Fitness) == min(Fitness)
        M = ones(Atom_Num,1) / Atom_Num;
    else
        M = exp(-(Fitness - max(Fitness)) ./ (max(Fitness) - min(Fitness)));
        M = M ./ sum(M);
    end

    G = exp(-20 * Iteration / Max_Iteration);

    Kbest = Atom_Num - (Atom_Num - 2) * (Iteration / Max_Iteration)^0.5;
    Kbest = floor(Kbest) + 1;

    [~, Index_M] = sort(M, 'descend');

    E = zeros(Atom_Num, Dim);

    MK = sum(Atom_Pop(Index_M(1:Kbest),:),1) / Kbest;

    for i = 1:Atom_Num

        Distance = norm(Atom_Pop(i,:) - MK);

        for k = 1:Kbest
            j = Index_M(k);
            Potential = LJPotential(Atom_Pop(i,:), Atom_Pop(j,:), ...
                Iteration, Max_Iteration, Distance);

            E(i,:) = E(i,:) + rand(1,Dim) * Potential .* ...
                ((Atom_Pop(j,:) - Atom_Pop(i,:)) / ...
                (norm(Atom_Pop(i,:) - Atom_Pop(j,:)) + eps));
        end

        E(i,:) = alpha * E(i,:) + beta * (X_Best - Atom_Pop(i,:));
    end

    a = E ./ (M + eps);
    Acc = a .* G;

end

%% ------------------------------------------------------------------------
function Potential = LJPotential(Atom1, Atom2, Iteration, Max_Iteration, s)

    r = norm(Atom1 - Atom2);
    c = (1 - (Iteration - 1) / Max_Iteration)^3;

    rsmin = 1.1 + 0.1 * sin(Iteration / Max_Iteration * pi / 2);
    rsmax = 1.24;

    if s == 0
        s = eps;
    end

    ratio = r / s;

    if ratio < rsmin
        rs = rsmin;
    elseif ratio > rsmax
        rs = rsmax;
    else
        rs = ratio;
    end

    Potential = c * (12 * (-rs)^(-13) - 6 * (-rs)^(-7));

end
