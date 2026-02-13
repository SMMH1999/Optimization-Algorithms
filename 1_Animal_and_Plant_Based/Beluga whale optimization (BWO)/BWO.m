function [bestFitness, bestPosition, convergenceCurve] = BWO(lb, ub, dim, nPop, maxItr, objFun)
    %==========================================================================
    % Beluga Whale Optimization (BWO)
    % Author: Original by [Author Unknown], Refactored for Benchmark Use
    %
    % Description:
    %   BWO is a nature-inspired metaheuristic algorithm based on the
    %   behavior of beluga whales. It balances exploration and exploitation
    %   using probabilistic updates and Levy flights.
    %
    % Inputs:
    %   lb       - Lower bound (scalar or 1xdim vector)
    %   ub       - Upper bound (scalar or 1xdim vector)
    %   dim      - Number of decision variables
    %   nPop     - Population size
    %   maxItr   - Maximum number of iterations
    %   objFun   - Handle to the objective function to minimize
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Position vector corresponding to bestFitness
    %   convergenceCurve - Array of bestFitness at each iteration
    %
    % Tunable Parameters:
    %   WF     - Probability of whale fall, decreases from 0.1 to 0.05
    %   KD     - Scale factor for Levy flight steps (0.05)
    %==========================================================================

    %% Initialization
    if isscalar(lb), lb = lb*ones(1, dim); end
    if isscalar(ub), ub = ub*ones(1, dim); end

    pos = rand(nPop, dim).*(ub-lb) + lb;   % Initial population
    fit = inf(nPop,1);
    Counts_run = 0;

    % Evaluate initial population
    for i = 1:nPop
        fit(i) = objFun(pos(i,:));
        Counts_run = Counts_run + 1;
    end

    [bestFitness, idx] = min(fit);
    bestPosition = pos(idx,:);

    convergenceCurve = inf(1, maxItr);

    %% Main Loop
    for t = 1:maxItr
        newPos = pos;

        WF = 0.1 - 0.05*(t/maxItr);                 % Whale fall probability
        kk = (1 - 0.5*t/maxItr)*rand(nPop,1);      % Exploration/Exploitation factor

        for i = 1:nPop
            if kk(i) > 0.5
                %% Exploration Phase
                r1 = rand(); r2 = rand();
                RJ = ceil(nPop*rand);
                while RJ == i
                    RJ = ceil(nPop*rand);
                end

                if dim <= nPop/5
                    params = randperm(dim,2);
                    newPos(i,params(1)) = pos(i,params(1)) + ...
                        (pos(RJ,params(1)) - pos(i,params(2)))*(r1+1)*sin(r2*360);
                    newPos(i,params(2)) = pos(i,params(2)) + ...
                        (pos(RJ,params(1)) - pos(i,params(2)))*(r1+1)*cos(r2*360);
                else
                    params = randperm(dim);
                    for j = 1:floor(dim/2)
                        newPos(i,2*j-1) = pos(i,params(2*j-1)) + ...
                            (pos(RJ,params(1)) - pos(i,params(2*j-1)))*(r1+1)*sin(r2*360);
                        newPos(i,2*j) = pos(i,params(2*j)) + ...
                            (pos(RJ,params(1)) - pos(i,params(2*j)))*(r1+1)*cos(r2*360);
                    end
                end

            else
                %% Exploitation Phase (Levy Flight)
                r3 = rand(); r4 = rand(); C1 = 2*r4*(1-t/maxItr);
                RJ = ceil(nPop*rand);
                while RJ == i
                    RJ = ceil(nPop*rand);
                end

                alpha = 3/2;
                sigma = (gamma(1+alpha)*sin(pi*alpha/2)/(gamma((1+alpha)/2)*alpha*2^((alpha-1)/2)))^(1/alpha);
                u = randn(1,dim).*sigma; v = randn(1,dim);
                LevyFlight = 0.05 * u ./ abs(v).^(1/alpha);

                newPos(i,:) = r3*bestPosition - r4*pos(i,:) + C1*LevyFlight.*(pos(RJ,:) - pos(i,:));
            end

            %% Boundary Handling
            newPos(i,:) = max(min(newPos(i,:), ub), lb);

            %% Fitness Update
            newFit = objFun(newPos(i,:));
            Counts_run = Counts_run + 1;
            if newFit < fit(i)
                pos(i,:) = newPos(i,:);
                fit(i) = newFit;
            end
        end

        %% Whale Fall Update
        for i = 1:nPop
            if kk(i) <= WF
                RJ = ceil(nPop*rand);
                while RJ == i
                    RJ = ceil(nPop*rand);
                end
                r5 = rand(); r6 = rand(); r7 = rand();
                C2 = 2*nPop*WF;
                stepSize = r7*(ub-lb).*exp(-C2*t/maxItr);
                newPos(i,:) = (r5*pos(i,:) - r6*pos(RJ,:)) + stepSize;

                %% Boundary Handling
                newPos(i,:) = max(min(newPos(i,:), ub), lb);

                %% Fitness Update
                newFit = objFun(newPos(i,:));
                Counts_run = Counts_run + 1;
                if newFit < fit(i)
                    pos(i,:) = newPos(i,:);
                    fit(i) = newFit;
                end
            end
        end

        %% Update Global Best
        [currBestFit, idx] = min(fit);
        if currBestFit < bestFitness
            bestFitness = currBestFit;
            bestPosition = pos(idx,:);
        end

        %% Store Convergence
        convergenceCurve(t) = bestFitness;
    end

end
