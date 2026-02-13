function [bestFitness, bestPosition, convergenceCurve] = FA(lb, ub, dim, nPop, maxItr, objFun, nonlcon)
    % =====================================================================
    % Firefly Algorithm (FA)
    % Author: Xin-She Yang, refactored by user
    % Description: Metaheuristic optimization algorithm inspired by the
    %              flashing behavior of fireflies.
    %
    % Inputs:
    %   lb        : Lower bounds (scalar or vector) [1 x dim]
    %   ub        : Upper bounds (scalar or vector) [1 x dim]
    %   dim       : Dimension of the problem
    %   nPop      : Number of fireflies (population size)
    %   maxItr    : Maximum number of iterations
    %   objFun    : Handle to objective function (minimization)
    %   nonlcon   : Optional handle to nonlinear constraints (return [g, geq])
    %
    % Outputs:
    %   bestFitness       : Best fitness value found
    %   bestPosition      : Position vector of the best solution
    %   convergenceCurve  : Vector of best fitness at each iteration
    %
    % Tunable parameters:
    %   alpha   : Randomness (0-1)
    %   gamma   : Light absorption coefficient
    %   betamin : Minimum attractiveness
    % =====================================================================

    if nargin < 7, nonlcon = []; end

    % ---- Parameters ----
    alpha = 0.25;        % Randomness
    betamin = 0.20;      % Minimum attractiveness
    gamma = 1.0;         % Light absorption coefficient

    % ---- Initialization ----
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end

    scale = abs(ub - lb);

    % Initialize fireflies
    ns = repmat(lb, nPop, 1) + rand(nPop, dim).*repmat(scale, nPop, 1);
    Lightn = inf(nPop,1);  % Fitness values
    convergenceCurve = zeros(maxItr,1);

    % ---- Main Loop ----
    for t = 1:maxItr
        % Evaluate fitness with optional nonlinear constraints
        for i = 1:nPop
            Lightn(i) = objFun(ns(i,:)) + penalty(nonlcon, ns(i,:));
        end

        % Rank fireflies (ascending, because minimization)
        [Lightn, idx] = sort(Lightn);
        ns = ns(idx,:);

        % Store the best solution
        bestPosition = ns(1,:);
        bestFitness = Lightn(1);
        convergenceCurve(t) = bestFitness;

        % Move fireflies
        ns = move_fireflies(ns, Lightn, alpha, betamin, gamma, lb, ub);

        % Reduce randomness
        alpha = alpha * 0.97;
    end

    % ------------------- Subfunctions -------------------

    function ns = move_fireflies(ns, Lightn, alpha, betamin, gamma, lb, ub)
        n = size(ns,1);
        dim = size(ns,2);
        ns_new = ns;

        for i = 1:n
            for j = 1:n
                if Lightn(i) > Lightn(j)
                    r = norm(ns(i,:) - ns(j,:));
                    beta = (1 - betamin)*exp(-gamma*r^2) + betamin;
                    ns_new(i,:) = ns_new(i,:)*(1-beta) + ns(j,:)*beta + alpha*(rand(1,dim)-0.5).*abs(ub-lb);
                end
            end
        end

        % Apply bounds
        ns_new = max(ns_new, lb);
        ns_new = min(ns_new, ub);
        ns = ns_new;

    end

    function z = penalty(nonlcon, x)
        z = 0;
        if isempty(nonlcon), return; end
        [g, geq] = nonlcon(x);
        lam = 1e15; lameq = 1e15;
        for k = 1:length(g), z = z + lam*g(k)^2*(g(k)>0); end
        for k = 1:length(geq), z = z + lameq*geq(k)^2*(geq(k)~=0); end
    end

end
