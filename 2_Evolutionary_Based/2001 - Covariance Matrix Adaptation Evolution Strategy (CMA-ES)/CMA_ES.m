function [bestFitness, bestPosition, convergenceCurve] = CMA_ES(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
    % Abbreviation: CMA-ES
    %
    % Original Developer:
    % S. Mostapha Kalami Heris (Yarpiz Team)
    %
    % Refactored to Benchmark Format for OAF Project
    %
    % Description:
    % Covariance Matrix Adaptation Evolution Strategy (CMA-ES) is a stochastic,
    % derivative-free optimization algorithm for solving nonlinear, nonconvex
    % continuous optimization problems. It adapts the covariance matrix of a
    % multivariate normal distribution to efficiently sample promising regions
    % of the search space.
    %
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables (dimension)
    %   nPop      - Population size (lambda: number of offspring)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Best decision vector found (1×dim)
    %   convergenceCurve - Best fitness value at each iteration (maxItr×1)
    %
    % Tunable Parameters (Derived Internally):
    %   mu       - Number of parents
    %   w        - Recombination weights
    %   mu_eff   - Effective selection mass
    %   cs, ds   - Step-size control parameters
    %   cc, c1, cmu - Covariance adaptation parameters
    %
    % Notes:
    %   - Minimization is assumed.
    %   - Boundary constraints are enforced via projection.
    %   - Original CMA-ES logic preserved exactly.
    % =========================================================================

    %% Bound Handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    %% CMA-ES Parameters

    lambda = nPop;                % Offspring size
    mu = round(lambda/2);         % Number of parents

    w = log(mu+0.5) - log(1:mu);
    w = w / sum(w);
    mu_eff = 1 / sum(w.^2);

    sigma0 = 0.3 * mean(ub - lb);

    cs = (mu_eff + 2) / (dim + mu_eff + 5);
    ds = 1 + cs + 2 * max(sqrt((mu_eff-1)/(dim+1)) - 1, 0);
    ENN = sqrt(dim) * (1 - 1/(4*dim) + 1/(21*dim^2));

    cc = (4 + mu_eff/dim) / (4 + dim + 2*mu_eff/dim);
    c1 = 2 / ((dim + 1.3)^2 + mu_eff);
    alpha_mu = 2;
    cmu = min(1 - c1, alpha_mu*(mu_eff-2+1/mu_eff)/((dim+2)^2 + alpha_mu*mu_eff/2));
    hth = (1.4 + 2/(dim+1)) * ENN;

    %% Initialization

    meanPos = lb + rand(1,dim) .* (ub - lb);
    meanStep = zeros(1,dim);

    ps = zeros(1,dim);
    pc = zeros(1,dim);
    C = eye(dim);
    sigma = sigma0;

    bestPosition = meanPos;
    bestFitness = objFun(bestPosition);

    convergenceCurve = zeros(maxItr,1);

    %% Main Loop

    for t = 1:maxItr

        % Generate Offspring
        popStep = zeros(lambda, dim);
        popPos  = zeros(lambda, dim);
        popFit  = zeros(lambda,1);

        for i = 1:lambda
            popStep(i,:) = mvnrnd(zeros(1,dim), C);
            popPos(i,:)  = meanPos + sigma * popStep(i,:);

            % Boundary control
            popPos(i,:) = max(popPos(i,:), lb);
            popPos(i,:) = min(popPos(i,:), ub);

            popFit(i) = objFun(popPos(i,:));

            if popFit(i) < bestFitness
                bestFitness = popFit(i);
                bestPosition = popPos(i,:);
            end
        end

        % Sort Population
        [popFit, idx] = sort(popFit);
        popStep = popStep(idx,:);
        popPos  = popPos(idx,:);

        convergenceCurve(t) = bestFitness;

        if t == maxItr
            break;
        end

        %% Update Mean
        newMeanStep = zeros(1,dim);
        for j = 1:mu
            newMeanStep = newMeanStep + w(j) * popStep(j,:);
        end

        newMeanPos = meanPos + sigma * newMeanStep;
        newMeanPos = max(newMeanPos, lb);
        newMeanPos = min(newMeanPos, ub);

        newMeanFit = objFun(newMeanPos);
        if newMeanFit < bestFitness
            bestFitness = newMeanFit;
            bestPosition = newMeanPos;
        end

        %% Step-Size Update
        ps = (1-cs)*ps + sqrt(cs*(2-cs)*mu_eff) * (newMeanStep / chol(C)');
        sigma = sigma * exp(cs/ds * (norm(ps)/ENN - 1))^0.3;

        %% Covariance Matrix Update
        if norm(ps)/sqrt(1-(1-cs)^(2*t)) < hth
            hs = 1;
        else
            hs = 0;
        end

        delta = (1-hs)*cc*(2-cc);
        pc = (1-cc)*pc + hs*sqrt(cc*(2-cc)*mu_eff)*newMeanStep;

        C = (1-c1-cmu)*C + c1*(pc'*pc + delta*C);

        for j = 1:mu
            C = C + cmu*w(j) * (popStep(j,:)' * popStep(j,:));
        end

        % Ensure Positive Semi-Definite
        [V,E] = eig(C);
        if any(diag(E) < 0)
            E = max(E,0);
            C = V * E / V;
        end

        meanPos = newMeanPos;
        meanStep = newMeanStep;

    end

end
