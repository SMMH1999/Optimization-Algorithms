function [bestFitness, bestPosition, convergenceCurve] = AO(lb, ub, dim, nPop, maxItr, objFun)
    % Aquila Optimizer (AO)
    %% Initialization

    % Ensure row vectors for bounds
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    % Initialize population
    X = rand(nPop, dim) .* (ub - lb) + lb;
    Xnew = zeros(nPop, dim);

    fitness = inf(nPop, 1);

    bestFitness = inf;
    bestPosition = zeros(1, dim);

    alpha = 0.1;
    delta = 0.1;

    convergenceCurve = zeros(1, maxItr);

    %% Main Loop
    for t = 1:maxItr

        % ---- Evaluation ----
        for i = 1:nPop

            % Boundary control
            X(i,:) = max(X(i,:), lb);
            X(i,:) = min(X(i,:), ub);

            % Fitness evaluation
            fitness(i) = objFun(X(i,:));

            % Update global best
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = X(i,:);
            end
        end

        % ---- Parameters ----
        G2 = 2 * rand() - 1;
        G1 = 2 * (1 - (t / maxItr));

        to = 1:dim;
        u = 0.0265;
        r0 = 10;
        r = r0 + u * to;
        omega = 0.005;
        phi0 = 3 * pi / 2;
        phi = -omega * to + phi0;
        x_spiral = r .* sin(phi);
        y_spiral = r .* cos(phi);

        QF = t^((2 * rand() - 1) / (1 - maxItr)^2);

        % ---- Position Update ----
        for i = 1:nPop

            if t <= (2/3) * maxItr

                if rand < 0.5
                    % Eq. (3)-(4)
                    Xnew(i,:) = bestPosition * (1 - t/maxItr) + ...
                        (mean(X(i,:)) - bestPosition) * rand;
                else
                    % Eq. (5)
                    levyStep = Levy(dim);
                    randomIndex = randi(nPop);
                    Xnew(i,:) = bestPosition .* levyStep + ...
                        X(randomIndex,:) + ...
                        (y_spiral - x_spiral) * rand;
                end

            else

                if rand < 0.5
                    % Eq. (13)
                    Xnew(i,:) = (bestPosition - mean(X)) * alpha ...
                        - rand + ...
                        ((ub - lb) * rand + lb) * delta;
                else
                    % Eq. (14)
                    levyStep = Levy(dim);
                    Xnew(i,:) = QF * bestPosition ...
                        - (G2 * X(i,:) * rand) ...
                        - G1 .* levyStep ...
                        + rand * G2;
                end

            end

            % Boundary control
            Xnew(i,:) = max(Xnew(i,:), lb);
            Xnew(i,:) = min(Xnew(i,:), ub);

            % Greedy selection
            newFitness = objFun(Xnew(i,:));
            if newFitness < fitness(i)
                X(i,:) = Xnew(i,:);
                fitness(i) = newFitness;

                if newFitness < bestFitness
                    bestFitness = newFitness;
                    bestPosition = Xnew(i,:);
                end
            end
        end

        convergenceCurve(t) = bestFitness;

    end

end

%% Levy Flight Function (Local)
function step = Levy(d)
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
        (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    step = u ./ abs(v).^(1/beta);
end
