function [bestFitness, bestPosition, convergenceCurve] = MFO(LB, UB, Dim, popSize, maxItr, Cost_Function)

    % Initialize moth positions
    mothPos = initialization(popSize, Dim, UB, LB);

    % Evaluate initial fitness
    mothFit = arrayfun(@(i) Cost_Function(mothPos(i,:)), 1:popSize).';

    % Convergence curve
    convergenceCurve = inf(maxItr, 1);

    % Initialize flames
    bestFlamesPos = mothPos;
    bestFlamesFit = mothFit;

    % Main loop
    for itr = 1:maxItr

        % Sort moths by fitness
        [mothFit, idx] = sort(mothFit);
        mothPos = mothPos(idx, :);

        % Number of flames (Eq. 3.14)
        flameNo = round(popSize - itr * ((popSize - 1) / maxItr));

        % Merge current moths with previous flames
        combinedPos = [mothPos; bestFlamesPos];
        combinedFit = [mothFit; bestFlamesFit];

        % Sort combined population
        [combinedFit, idxCombined] = sort(combinedFit);
        combinedPos = combinedPos(idxCombined, :);

        % Update flames
        bestFlamesPos = combinedPos(1:popSize, :);
        bestFlamesFit = combinedFit(1:popSize);

        % Best solution so far
        bestFitness  = bestFlamesFit(1);
        bestPosition = bestFlamesPos(1, :);

        % Decreasing parameter a
        a = -1 - itr / maxItr;

        % Update moth positions
        for i = 1:popSize

            flameIdx = min(i, flameNo);

            for d = 1:Dim
                dist = abs(bestFlamesPos(flameIdx, d) - mothPos(i, d));
                t = (a - 1) * rand + 1;

                mothPos(i, d) = dist * exp(t) * cos(2*pi*t) ...
                              + bestFlamesPos(flameIdx, d);
            end
        end

        % Boundary control
        mothPos = max(min(mothPos, UB), LB);

        % Recalculate fitness
        mothFit = arrayfun(@(i) Cost_Function(mothPos(i,:)), 1:popSize).';

        % Update convergence curve
        if itr == 1
            convergenceCurve(itr) = bestFitness;
        else
            convergenceCurve(itr) = min(convergenceCurve(itr-1), bestFitness);
        end
    end

    % Safety: enforce monotonic convergence
    convergenceCurve = cummin(convergenceCurve);
end

%% ====================================================================
% Generate random population inside bounds
function X = initialization(nAgents, dim, ub, lb)

    if isscalar(ub) && isscalar(lb)
        X = rand(nAgents, dim) .* (ub - lb) + lb;
    else
        X = zeros(nAgents, dim);
        for d = 1:dim
            X(:, d) = rand(nAgents, 1) .* (ub(d) - lb(d)) + lb(d);
        end
    end
end
