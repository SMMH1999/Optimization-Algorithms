function [bestFitness, bestPosition, convergenceCurve] = FPA(lb, ub, dim, nPop, maxItr, objFun)
    % -------------------------------------------------------------------------
    % Flower Pollination Algorithm (FPA)
    % Author: Xin-She Yang (2012)
    %
    % Description:
    %   A nature-inspired global optimization algorithm based on the
    %   pollination process of flowers. Supports single-objective problems.
    %
    % Inputs:
    %   lb        - lower bounds (1xd vector or scalar)
    %   ub        - upper bounds (1xd vector or scalar)
    %   dim       - number of decision variables (dimension)
    %   nPop      - population size (number of flowers/pollens)
    %   maxItr    - maximum number of iterations
    %   objFun    - handle to the objective function (fitness function)
    %
    % Outputs:
    %   bestFitness       - best objective function value found
    %   bestPosition      - position (solution vector) corresponding to bestFitness
    %   convergenceCurve  - bestFitness value at each iteration
    %
    % Tunable Parameters:
    %   p = 0.8   - probability switch between global and local pollination
    %
    % -------------------------------------------------------------------------

    % Probability switch for global/local pollination
    p = 0.8;

    % Ensure lb and ub are vectors of correct dimension
    if isscalar(lb), lb = lb*ones(1, dim); end
    if isscalar(ub), ub = ub*ones(1, dim); end

    % Initialize population and fitness
    Sol = repmat(lb, nPop, 1) + rand(nPop, dim) .* repmat((ub - lb), nPop, 1);
    Fitness = objFun(Sol);

    % Find the initial best
    [bestFitness, idx] = min(Fitness);
    bestPosition = Sol(idx,:);
    convergenceCurve = zeros(1, maxItr);

    % Main loop
    for t = 1:maxItr
        for i = 1:nPop
            if rand > p
                % Global pollination (Levy flights)
                L = Levy(dim);
                S = Sol(i,:) + L .* (Sol(i,:) - bestPosition);
            else
                % Local pollination
                epsilon = rand;
                JK = randperm(nPop,2);
                S = Sol(i,:) + epsilon * (Sol(JK(1),:) - Sol(JK(2),:));
            end

            % Enforce bounds
            S = simplebounds(S, lb, ub);

            % Evaluate new solution
            Fnew = objFun(S);

            % Accept new solution if better
            if Fnew <= Fitness(i)
                Sol(i,:) = S;
                Fitness(i) = Fnew;
            end

            % Update global best
            if Fnew <= bestFitness
                bestPosition = S;
                bestFitness = Fnew;
            end
        end
        % Record convergence
        convergenceCurve(t) = bestFitness;
    end

end

%% ---------------------- Local Helper Functions ------------------------

function s = simplebounds(s, lb, ub)
    % Apply lower and upper bounds
    s = max(s, lb);
    s = min(s, ub);
end

function L = Levy(d)
    % Generate Levy flight steps
    beta = 3/2;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
        (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1,d) * sigma;
    v = randn(1,d);
    step = u ./ abs(v).^(1/beta);
    L = 0.01 * step;
end
