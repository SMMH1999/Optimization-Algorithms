function [bestFitness, bestPosition, convergenceCurve] = CSO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Competitive Swarm Optimization (CSO)
    % =========================================================================
    % Author: Original CSO concept adapted and unified for benchmark format
    %
    % Description:
    % Competitive Swarm Optimization (CSO) is a population-based metaheuristic
    % algorithm inspired by swarm intelligence. Each particle updates its
    % velocity based on its personal best position and the global best
    % position with predefined probabilistic learning rates.
    %
    % =========================================================================
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Problem dimension (integer)
    %   nPop      - Population size (number of particles)
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Best solution vector found (1×dim)
    %   convergenceCurve - Best fitness at each iteration (maxItr×1)
    %
    % =========================================================================
    % Tunable Parameters (fixed internally to preserve original behavior):
    %   pbest_rate  = 0.5   % Probability of learning from personal best
    %   gbest_rate  = 0.5   % Probability of learning from global best
    %   alpha       = 0.5   % Learning coefficient
    %
    % =========================================================================
    % Structure:
    %   1) Initialization
    %   2) Fitness Evaluation
    %   3) Main Optimization Loop
    %   4) Personal & Global Best Update
    %   5) Convergence Curve Update
    % =========================================================================

    %% -------------------- Parameter Definition ----------------------------
    pbest_rate = 0.5;
    gbest_rate = 0.5;
    alpha      = 0.5;

    %% -------------------- Boundary Handling -------------------------------
    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end

    %% -------------------- Initialization ----------------------------------
    swarm    = lb + (ub - lb) .* rand(nPop, dim);
    velocity = zeros(nPop, dim);

    pbest    = swarm;
    pbestFit = zeros(nPop,1);

    for i = 1:nPop
        pbestFit(i) = objFun(swarm(i,:));
    end

    [bestFitness, idx] = min(pbestFit);
    bestPosition       = swarm(idx,:);

    convergenceCurve = zeros(maxItr,1);

    %% -------------------- Main Optimization Loop --------------------------
    for t = 1:maxItr

        % -------- Velocity Update --------
        for i = 1:nPop
            if rand < pbest_rate
                velocity(i,:) = velocity(i,:) + alpha * (pbest(i,:) - swarm(i,:));
            end
            if rand < gbest_rate
                velocity(i,:) = velocity(i,:) + alpha * (bestPosition - swarm(i,:));
            end
        end

        % -------- Position Update --------
        swarm = swarm + velocity;

        % -------- Boundary Control --------
        swarm = max(swarm, lb);
        swarm = min(swarm, ub);

        % -------- Fitness Evaluation --------
        fitness = zeros(nPop,1);
        for i = 1:nPop
            fitness(i) = objFun(swarm(i,:));
        end

        % -------- Personal Best Update --------
        for i = 1:nPop
            if fitness(i) < pbestFit(i)
                pbest(i,:)   = swarm(i,:);
                pbestFit(i)  = fitness(i);
            end
        end

        % -------- Global Best Update --------
        [currentBestFit, idx] = min(pbestFit);
        if currentBestFit < bestFitness
            bestFitness  = currentBestFit;
            bestPosition = pbest(idx,:);
        end

        % -------- Convergence Curve --------
        convergenceCurve(t) = bestFitness;
    end

end
