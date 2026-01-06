function [bestScore, bestPos, curve] = GWO(LB, UB, Dim, populationNo, maxItr, objective)
    % GWO  Grey Wolf Optimizer (minimization)
    %
    % Inputs:
    %   LB, UB         - lower/upper bounds (scalar or 1xDim)
    %   Dim            - problem dimension
    %   populationNo   - number of search agents
    %   maxItr         - maximum iterations
    %   objective      - function handle, objective(x) -> scalar
    %
    % Outputs:
    %   bestScore      - best objective value found
    %   bestPos        - best position (1xDim)
    %   curve          - convergence curve (maxItr x 1)

    % Normalize bounds
    if isscalar(LB), LB = repmat(LB, 1, Dim); end
    if isscalar(UB), UB = repmat(UB, 1, Dim); end

    % Leaders
    alphaPos  = zeros(1, Dim);
    betaPos   = zeros(1, Dim);
    deltaPos  = zeros(1, Dim);

    alphaScore = inf;
    betaScore  = inf;
    deltaScore = inf;

    % Initialize population
    Positions = rand(populationNo, Dim) .* (UB - LB) + LB;

    curve = zeros(maxItr, 1);

    for it = 1:maxItr

        % Evaluate population & update leaders
        for i = 1:populationNo
            % Keep inside bounds
            Positions(i, :) = min(max(Positions(i, :), LB), UB);

            fitness = objective(Positions(i, :));

            if fitness < alphaScore
                alphaScore = fitness;
                alphaPos   = Positions(i, :);
            elseif fitness < betaScore
                betaScore = fitness;
                betaPos   = Positions(i, :);
            elseif fitness < deltaScore
                deltaScore = fitness;
                deltaPos   = Positions(i, :);
            end
        end

        % a decreases linearly from 2 to 0
        a = 2 - it * (2 / maxItr);

        % Update positions
        for i = 1:populationNo
            for j = 1:Dim
                A1 = 2 * a * rand() - a;
                C1 = 2 * rand();
                D_alpha = abs(C1 * alphaPos(j) - Positions(i, j));
                X1 = alphaPos(j) - A1 * D_alpha;

                A2 = 2 * a * rand() - a;
                C2 = 2 * rand();
                D_beta = abs(C2 * betaPos(j) - Positions(i, j));
                X2 = betaPos(j) - A2 * D_beta;

                A3 = 2 * a * rand() - a;
                C3 = 2 * rand();
                D_delta = abs(C3 * deltaPos(j) - Positions(i, j));
                X3 = deltaPos(j) - A3 * D_delta;

                Positions(i, j) = (X1 + X2 + X3) / 3;
            end
        end

        curve(it) = alphaScore;
    end

    bestScore = alphaScore;
    bestPos   = alphaPos;
end
