function [bestFitness, bestPosition, convergenceCurve] = COA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Coyote Optimization Algorithm (COA)
    % =========================================================================
    % Author: Juliano Pierezan and Leandro dos Santos Coelho (2018)
    % Refactored to benchmark single-function format.
    %
    % Reference:
    % Pierezan, J. and Coelho, L. S. (2018). "Coyote Optimization Algorithm:
    % A new metaheuristic for global optimization problems",
    % IEEE Congress on Evolutionary Computation (CEC), pp. 2633–2640.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        - Lower bound (1×dim or scalar)
    %   ub        - Upper bound (1×dim or scalar)
    %   dim       - Number of decision variables
    %   nPop      - Total population size (n_packs × n_coyotes)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective function value found
    %   bestPosition     - Decision variables of the best solution
    %   convergenceCurve - Best fitness value at each iteration
    %
    % Tunable Parameters (internal structure):
    %   n_packs   - Number of packs (fixed = 5)
    %   n_coy     - Coyotes per pack = nPop / n_packs
    %   p_leave   - Probability of leaving a pack
    %   Ps        - Probability parameter for reproduction
    %
    % Notes:
    %   - Minimization problem
    %   - Boundary constraints strictly enforced
    %   - Fully benchmark-compatible structure
    % =========================================================================

    %% ------------------------- Initialization -------------------------------

    if numel(lb) == 1
        VarMin = lb * ones(1, dim);
    else
        VarMin = lb;
    end

    if numel(ub) == 1
        VarMax = ub * ones(1, dim);
    else
        VarMax = ub;
    end

    % Pack structure (fixed to preserve original logic)
    n_packs = 5;
    n_coy   = floor(nPop / n_packs);
    nPop    = n_packs * n_coy;

    p_leave = 0.005 * n_coy^2;
    Ps      = 1 / dim;

    % Population initialization (Eq. 2)
    coyotes = repmat(VarMin, nPop, 1) + ...
        rand(nPop, dim) .* (repmat(VarMax - VarMin, nPop, 1));

    ages    = zeros(nPop, 1);
    packs   = reshape(randperm(nPop), n_packs, []);
    costs   = zeros(nPop, 1);

    %% ----------------------- Fitness Evaluation -----------------------------

    for i = 1:nPop
        costs(i) = objFun(coyotes(i, :));
    end

    [bestFitness, idx] = min(costs);
    bestPosition = coyotes(idx, :);

    convergenceCurve = zeros(maxItr, 1);

    %% --------------------------- Main Loop ----------------------------------

    for t = 1:maxItr

        for p = 1:n_packs

            packIdx     = packs(p, :);
            coy_aux     = coyotes(packIdx, :);
            cost_aux    = costs(packIdx);
            age_aux     = ages(packIdx);

            % Sort coyotes (Eq. 5)
            [cost_aux, order] = sort(cost_aux);
            coy_aux = coy_aux(order, :);
            age_aux = age_aux(order);

            c_alpha  = coy_aux(1, :);
            tendency = median(coy_aux, 1); % Eq. 6

            % Social update (Eq. 12–14)
            for c = 1:n_coy

                rc1 = c;
                while rc1 == c
                    rc1 = randi(n_coy);
                end

                rc2 = c;
                while rc2 == c || rc2 == rc1
                    rc2 = randi(n_coy);
                end

                new_c = coy_aux(c, :) + ...
                    rand * (c_alpha - coy_aux(rc1, :)) + ...
                    rand * (tendency - coy_aux(rc2, :));

                new_c = min(max(new_c, VarMin), VarMax);

                new_cost = objFun(new_c);

                if new_cost < cost_aux(c)
                    coy_aux(c, :) = new_c;
                    cost_aux(c)   = new_cost;
                end
            end

            %% Reproduction (Eq. 7 + Alg. 1)
            parents = randperm(n_coy, 2);

            prob1 = (1 - Ps) / 2;
            prob2 = prob1;

            pdr = randperm(dim);
            p1  = zeros(1, dim);
            p2  = zeros(1, dim);

            p1(pdr(1)) = 1;
            p2(pdr(2)) = 1;

            r = rand(1, dim - 2);
            p1(pdr(3:end)) = r < prob1;
            p2(pdr(3:end)) = r > 1 - prob2;

            n_mask = ~(p1 | p2);

            pup = p1 .* coy_aux(parents(1), :) + ...
                p2 .* coy_aux(parents(2), :) + ...
                n_mask .* (VarMin + rand(1, dim) .* (VarMax - VarMin));

            pup = min(max(pup, VarMin), VarMax);
            pup_cost = objFun(pup);

            worse = find(pup_cost < cost_aux);

            if ~isempty(worse)
                [~, olderIdx] = max(age_aux(worse));
                replaceIdx = worse(olderIdx);

                coy_aux(replaceIdx, :) = pup;
                cost_aux(replaceIdx)   = pup_cost;
                age_aux(replaceIdx)    = 0;
            end

            % Update pack
            coyotes(packIdx, :) = coy_aux;
            costs(packIdx)      = cost_aux;
            ages(packIdx)       = age_aux;
        end

        %% Pack exchange (Eq. 4)
        if n_packs > 1
            if rand < p_leave
                rp = randperm(n_packs, 2);
                rc = randi(n_coy, 1, 2);

                tmp = packs(rp(1), rc(1));
                packs(rp(1), rc(1)) = packs(rp(2), rc(2));
                packs(rp(2), rc(2)) = tmp;
            end
        end

        %% Age update
        ages = ages + 1;

        %% Global best update
        [currentBest, idx] = min(costs);
        if currentBest < bestFitness
            bestFitness  = currentBest;
            bestPosition = coyotes(idx, :);
        end

        convergenceCurve(t) = bestFitness;
    end

end
