function [bestFitness, bestPosition, convergenceCurve] = MRFO(lb, ub, dim, nPop, maxItr, objFun)
    %--------------------------------------------------------------------------
    % MRFO - Manta Ray Foraging Optimization
    %--------------------------------------------------------------------------
    % Developed by: Original code by W. Zhao, Z. Zhang, L. Wang (2019)
    % Refactored for MATLAB benchmark usage
    %
    % Description:
    %   MRFO is a bio-inspired optimization algorithm that mimics the foraging
    %   behavior of manta rays. It uses chain foraging, cyclone foraging, and
    %   somersault foraging mechanisms to explore and exploit the search space.
    %
    % Inputs:
    %   lb       - Lower bound of search space (scalar or 1xDim vector)
    %   ub       - Upper bound of search space (scalar or 1xDim vector)
    %   dim      - Number of decision variables (problem dimensionality)
    %   nPop     - Population size (number of manta rays)
    %   maxItr   - Maximum number of iterations
    %   objFun   - Handle to objective function: fitness = objFun(position)
    %
    % Outputs:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Best solution vector found
    %   convergenceCurve - Vector of best fitness at each iteration
    %
    % Tunable Parameters:
    %   - Chain foraging weight: implicitly handled via algorithm updates
    %   - Cyclone foraging weight: implicitly handled via algorithm updates
    %   - Somersault factor (S): set to 2 in the original algorithm
    %--------------------------------------------------------------------------

    %% Initialization
    if isscalar(lb), lb = lb*ones(1, dim); end
    if isscalar(ub), ub = ub*ones(1, dim); end

    PopPos = rand(nPop, dim).*(ub - lb) + lb;            % Initialize population
    PopFit = arrayfun(@(i) objFun(PopPos(i,:)), 1:nPop); % Evaluate fitness

    [bestFitness, idx] = min(PopFit);                   % Best initialization
    bestPosition = PopPos(idx,:);

    convergenceCurve = zeros(maxItr,1);                 % History of best fitness

    S = 2;  % Somersault factor

    %% Main Loop
    for t = 1:maxItr
        coef = t / maxItr;
        newPopPos = zeros(nPop, dim);

        % Update first individual
        if rand < 0.5
            r1 = rand;
            Beta = 2*exp(r1*((maxItr-t+1)/maxItr)) * sin(2*pi*r1);
            if coef > rand
                newPopPos(1,:) = bestPosition + rand(1,dim).*(bestPosition - PopPos(1,:)) + Beta*(bestPosition - PopPos(1,:));
            else
                indivRand = rand(1,dim).*(ub-lb) + lb;
                newPopPos(1,:) = indivRand + rand(1,dim).*(indivRand - PopPos(1,:)) + Beta*(indivRand - PopPos(1,:));
            end
        else
            Alpha = 2*rand(1,dim).*(-log(rand(1,dim))).^0.5;
            newPopPos(1,:) = PopPos(1,:) + rand(1,dim).*(bestPosition - PopPos(1,:)) + Alpha.*(bestPosition - PopPos(1,:));
        end

        % Update remaining individuals
        for i = 2:nPop
            if rand < 0.5
                r1 = rand;
                Beta = 2*exp(r1*((maxItr-t+1)/maxItr)) * sin(2*pi*r1);
                if coef > rand
                    newPopPos(i,:) = bestPosition + rand(1,dim).*(PopPos(i-1,:) - PopPos(i,:)) + Beta*(bestPosition - PopPos(i,:));
                else
                    indivRand = rand(1,dim).*(ub-lb) + lb;
                    newPopPos(i,:) = indivRand + rand(1,dim).*(PopPos(i-1,:) - PopPos(i,:)) + Beta*(indivRand - PopPos(i,:));
                end
            else
                Alpha = 2*rand(1,dim).*(-log(rand(1,dim))).^0.5;
                newPopPos(i,:) = PopPos(i,:) + rand(1,dim).*(PopPos(i-1,:) - PopPos(i,:)) + Alpha.*(bestPosition - PopPos(i,:));
            end
        end

        % Boundary check and fitness update
        for i = 1:nPop
            newPopPos(i,:) = SpaceBound(newPopPos(i,:), lb, ub);
            newFit = objFun(newPopPos(i,:));
            if newFit < PopFit(i)
                PopFit(i) = newFit;
                PopPos(i,:) = newPopPos(i,:);
            end
        end

        % Somersault foraging
        for i = 1:nPop
            newPopPos(i,:) = PopPos(i,:) + S*(rand*bestPosition - rand*PopPos(i,:));
            newPopPos(i,:) = SpaceBound(newPopPos(i,:), lb, ub);
            newFit = objFun(newPopPos(i,:));
            if newFit < PopFit(i)
                PopFit(i) = newFit;
                PopPos(i,:) = newPopPos(i,:);
            end
        end

        % Update global best
        [currentBest, idx] = min(PopFit);
        if currentBest < bestFitness
            bestFitness = currentBest;
            bestPosition = PopPos(idx,:);
        end

        convergenceCurve(t) = bestFitness;
    end

end

%% Local function: boundary handling
function X = SpaceBound(X, lb, ub)
    X = min(max(X, lb), ub);
end
