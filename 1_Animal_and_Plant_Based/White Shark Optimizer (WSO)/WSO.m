function [bestFitness, bestPosition, convergenceCurve] = WSO(lb, ub, dim, nPop, maxItr, objFun)
    %==========================================================================
    % White Shark Optimizer (WSO) - Refactored Benchmark Version
    %
    % Author: Malik Braik
    % Reference: Malik Braik et al., "White Shark Optimizer: A novel bio-inspired
    % meta-heuristic algorithm for global optimization problems", Knowledge-Based Systems, 2022.
    %
    %--------------------------------------------------------------------------
    % INPUTS:
    % lb       : Lower bound (scalar or 1xdim vector)
    % ub       : Upper bound (scalar or 1xdim vector)
    % dim      : Problem dimension (integer)
    % nPop     : Number of white sharks (population size)
    % maxItr   : Maximum number of iterations
    % objFun   : Handle to objective function
    %
    % OUTPUTS:
    % bestFitness      : Best objective function value found
    % bestPosition     : Best position vector found
    % convergenceCurve : Best fitness value at each iteration
    %
    %--------------------------------------------------------------------------
    % Tunable Parameters:
    % fmax, fmin  : Max/min frequency of wavy motion
    % tau         : Parameter for velocity update
    % a0, a1, a2  : Control parameters for movement
    % pmin, pmax  : Influence parameters for velocity
    %==========================================================================

    %% Initialize
    convergenceCurve = zeros(1,maxItr);
    WSO_Positions = initializePopulation(nPop, dim, ub, lb);  % Initial positions
    v = zeros(nPop, dim);                                     % Initial velocity
    fitness = zeros(nPop,1);

    % Evaluate initial population
    for i = 1:nPop
        fitness(i) = objFun(WSO_Positions(i,:));
    end

    [bestFitness, idx] = min(fitness);
    wbest = WSO_Positions;           % Individual best positions
    bestPosition = WSO_Positions(idx,:); % Global best

    %% WSO Parameters
    fmax = 0.75; fmin = 0.07;
    tau  = 4.11;
    mu   = 2 / abs(2 - tau - sqrt(tau^2 - 4*tau));
    pmin = 0.5; pmax = 1.5;
    a0   = 6.250; a1 = 100; a2 = 0.0005;

    %% Main Loop
    for ite = 1:maxItr
        mv  = 1 / (a0 + exp((maxItr/2 - ite)/a1));
        s_s = abs(1 - exp(-a2*ite/maxItr));

        p1 = pmax + (pmax-pmin)*exp(-(4*ite/maxItr)^2);
        p2 = pmin + (pmax-pmin)*exp(-(4*ite/maxItr)^2);

        %% Update velocities
        nu = randi(nPop, 1, nPop);
        for i = 1:nPop
            rr = 1 + rand()*2;  % rmin=1, rmax=3
            wr = abs((2*rand - (rand+rand))/rr);
            v(i,:) = mu*v(i,:) + wr*(wbest(nu(i),:) - WSO_Positions(i,:));
        end

        %% Update positions based on wavy motion
        for i = 1:nPop
            f = bestFitness + (fmax - fmin)/(fmax + fmin);
            a = WSO_Positions(i,:) > ub;
            b = WSO_Positions(i,:) < lb;
            wo = xor(a,b);

            if rand < mv
                WSO_Positions(i,:) = WSO_Positions(i,:).*(~wo) + (ub.*a + lb.*b);
            else
                WSO_Positions(i,:) = WSO_Positions(i,:) + v(i,:)/f;
            end
        end

        %% Fishing school influence
        for i = 1:nPop
            for j = 1:dim
                if rand < s_s
                    Dist = abs(rand*(bestPosition(j) - WSO_Positions(i,j)));
                    if i == 1
                        WSO_Positions(i,j) = bestPosition(j) + rand*Dist*sign(rand-0.5);
                    else
                        WSO_Positions(i,j) = bestPosition(j) + rand*Dist*sign(rand-0.5);
                        WSO_Positions(i,j) = (WSO_Positions(i,j) + WSO_Positions(i-1,j))/2*rand;
                    end
                end
            end
        end

        %% Evaluate and update best positions
        for i = 1:nPop
            % Enforce bounds
            WSO_Positions(i,:) = max(min(WSO_Positions(i,:), ub), lb);

            fitVal = objFun(WSO_Positions(i,:));
            if fitVal < fitness(i)
                fitness(i) = fitVal;
                wbest(i,:) = WSO_Positions(i,:);
            end

            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = wbest(i,:);
            end
        end

        %% Store convergence
        convergenceCurve(ite) = bestFitness;
    end
end

%% Local function: Initialize Population
function pos = initializePopulation(nPop, dim, ub, lb)
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end
    pos = rand(nPop, dim).*(ub - lb) + lb;
end
