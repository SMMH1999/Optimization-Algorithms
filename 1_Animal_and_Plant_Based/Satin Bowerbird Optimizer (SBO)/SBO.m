function [bestFitness, bestPosition, convergenceCurve] = SBO(LB, UB, Dim, popSize, maxItr, Cost_Function)

    %% Parameters
    alpha     = 0.94;    % greatest step size
    pMutation = 0.05;    % mutation probability
    Z         = 0.02;    % sigma ratio

    % Expand bounds if scalar
    if isscalar(UB)
        UB = repmat(UB, 1, Dim);
        LB = repmat(LB, 1, Dim);
    end

    sigma = Z * (max(UB) - min(LB));   % Gaussian sigma

    %% Initialization
    X   = popgen(popSize, Dim, LB, UB);
    fit = arrayfun(@(i) Cost_Function(X(i,:)), 1:popSize);

    [bestFitness, idx] = min(fit);
    bestPosition = X(idx, :);

    convergenceCurve = zeros(1, maxItr);
    convergenceCurve(1) = bestFitness;

    %% Main loop
    for t = 2:maxItr

        % Probability weights (better fitness -> higher probability)
        F = 1 ./ (1 + max(0, fit));
        P = F ./ sum(F);

        % Generate new population
        Xnew   = X;
        fitNew = zeros(1, popSize);

        for i = 1:popSize
            j = rouletteSelect(P);   % selected guide once per agent

            lambda = alpha / (1 + P(j));

            for d = 1:Dim
                Xnew(i,d) = X(i,d) ...
                    + lambda * ( ((X(j,d) + bestPosition(d)) / 2) - X(i,d) );

                % Mutation
                if rand < pMutation
                    Xnew(i,d) = Xnew(i,d) + sigma * randn();
                end
            end
        end

        % Boundary control
        Xnew = max(min(Xnew, UB), LB);

        % Evaluate new solutions
        fitNew = arrayfun(@(i) Cost_Function(Xnew(i,:)), 1:popSize);

        % Combine and select elites
        X   = [X; Xnew];
        fit = [fit, fitNew];

        [fit, order] = sort(fit);
        X   = X(order(1:popSize), :);
        fit = fit(1:popSize);

        % Update global best
        if fit(1) < bestFitness
            bestFitness  = fit(1);
            bestPosition = X(1,:);
        end

        % Update convergence curve
        convergenceCurve(t) = min(convergenceCurve(t-1), bestFitness);
    end
end

%====================================================================
% Roulette wheel selection
function idx = rouletteSelect(P)
    C = cumsum(P);
    idx = find(rand <= C, 1, 'first');
end

%====================================================================
% Population generator
function X = popgen(n, d, LB, UB)
    X = LB + rand(n, d) .* (UB - LB);
end
