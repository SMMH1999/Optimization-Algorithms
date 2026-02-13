function [bestFitness, bestPosition, convergenceCurve] = AO(lb, ub, dim, nPop, maxIter, objFun)
    %==================================================================
    % 🌿 Artemisinin Optimization (AO) 🌿
    % Author(s): Chong Yuan, Dong Zhao, Ali Asghar Heidari, Lei Liu, Yi Chen, Zongda Wu, Huiling Chen
    % Reference: "Artemisinin Optimization based on Malaria Therapy:
    % Algorithm and Applications to Medical Image Segmentation", Displays, Elsevier, 2024
    %
    % Description:
    %   AO is a bio-inspired optimization algorithm based on malaria therapy
    %   mechanisms. This implementation supports minimization problems.
    %
    % Inputs:
    %   lb       : lower bound (scalar or 1 x dim vector)
    %   ub       : upper bound (scalar or 1 x dim vector)
    %   dim      : number of decision variables
    %   nPop     : population size
    %   maxIter  : maximum number of iterations
    %   objFun   : handle to objective function, e.g., @(x) sum(x.^2)
    %
    % Outputs:
    %   bestFitness      : best fitness value found
    %   bestPosition     : position of best solution
    %   convergenceCurve : best fitness value at each iteration
    %
    %==================================================================

    %% Initialization
    pop = initializePopulation(nPop, dim, ub, lb);  % Initial population

    % Evaluate initial fitness
    Fitness = zeros(1, nPop);
    for i = 1:nPop
        Fitness(i) = objFun(pop(i,:));
    end

    % Initial best solution
    [bestFitness, idx] = min(Fitness);
    bestPosition = pop(idx,:);

    % Containers
    New_pop = zeros(nPop, dim);
    Fitnorm = zeros(1, nPop);
    convergenceCurve = zeros(1, maxIter);

    %% Main loop
    for iter = 1:maxIter

        K = 1 - (iter^(1/6) / maxIter^(1/6));
        E = exp(-4*(iter/maxIter));

        % Update each agent
        for i = 1:nPop
            Fitnorm(i) = (Fitness(i) - min(Fitness)) / (max(Fitness) - min(Fitness) + eps);

            for j = 1:dim
                if rand < K
                    if rand < 0.5
                        New_pop(i,j) = pop(i,j) + E * pop(i,j) * (-1)^iter;
                    else
                        New_pop(i,j) = pop(i,j) + E * bestPosition(j) * (-1)^iter;
                    end
                else
                    New_pop(i,j) = pop(i,j);
                end

                if rand < Fitnorm(i)
                    A = randperm(nPop);
                    beta = rand/2 + 0.1;
                    New_pop(i,j) = pop(A(3),j) + beta * (pop(A(1),j) - pop(A(2),j));
                end
            end

            % Apply mutation and boundary reset
            New_pop(i,:) = mutationOperator(New_pop(i,:), pop(i,:), bestPosition, dim);
            New_pop(i,:) = boundaryReset(New_pop(i,:), ub, lb, dim, bestPosition);

            % Evaluate new solution
            tFitness = objFun(New_pop(i,:));

            % Update if better
            if tFitness < Fitness(i)
                pop(i,:) = New_pop(i,:);
                Fitness(i) = tFitness;
            end
        end

        % Update global best
        [fmin, idx] = min(Fitness);
        if fmin < bestFitness
            bestPosition = pop(idx,:);
            bestFitness = fmin;
        end

        % Record convergence
        convergenceCurve(iter) = bestFitness;
    end

end

%% =================== Helper Functions ===================

function Positions = initializePopulation(nAgents, dim, ub, lb)
    % Initializes the population within bounds
    if isscalar(ub)
        Positions = rand(nAgents, dim) * (ub - lb) + lb;
    else
        Positions = zeros(nAgents, dim);
        for i = 1:dim
            Positions(:,i) = rand(nAgents,1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end

function z = mutationOperator(z, x, best, dim)
    % Mutation based on probabilities
    for j = 1:dim
        if rand < 0.05
            z(j) = x(j);
        end
        if rand < 0.2
            z(j) = best(j);
        end
    end
end

function z = boundaryReset(z, ub, lb, dim, best)
    % Reset positions exceeding boundaries to best solution
    for j = 1:dim
        if z(j) > ub || z(j) < lb
            z(j) = best(j);
        end
    end
end
