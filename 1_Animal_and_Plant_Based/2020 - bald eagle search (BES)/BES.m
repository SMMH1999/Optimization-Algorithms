function [bestFitness, bestPosition, convergenceCurve] = BES(lb, ub, dim, nPop, maxItr, objFun)
    % Bald Eagle Search Optimization Algorithm (BES)
    % Author: Alsattar, H. A., Zaidan, A. A., & Zaidan, B. B. (2020)
    %
    % Description:
    %   BES is a meta-heuristic algorithm inspired by bald eagle hunting behavior.
    %   It iteratively searches and exploits the solution space to find the minimum
    %   of a given objective function.
    %
    % Inputs:
    %   lb      - scalar or 1xdim vector, lower bounds of decision variables
    %   ub      - scalar or 1xdim vector, upper bounds of decision variables
    %   dim     - integer, number of decision variables
    %   nPop    - integer, population size
    %   maxItr  - integer, maximum number of iterations
    %   objFun  - function handle, objective function to minimize
    %
    % Outputs:
    %   bestFitness       - best objective function value found
    %   bestPosition      - position vector corresponding to bestFitness
    %   convergenceCurve  - vector of bestFitness at each iteration

    % Initialize population and best solution
    pop.pos = lb + (ub-lb).*rand(nPop, dim);
    pop.cost = zeros(nPop,1);
    bestFitness = inf;
    bestPosition = zeros(1,dim);

    for i = 1:nPop
        pop.cost(i) = objFun(pop.pos(i,:));
        if pop.cost(i) < bestFitness
            bestFitness = pop.cost(i);
            bestPosition = pop.pos(i,:);
        end
    end

    convergenceCurve = zeros(1, maxItr);

    % Main loop
    for t = 1:maxItr
        [pop, bestFitness, bestPosition] = select_space(pop, nPop, bestFitness, bestPosition, lb, ub, dim, objFun);
        [pop, bestFitness, bestPosition] = search_space(pop, nPop, bestFitness, bestPosition, lb, ub, dim, objFun);
        [pop, bestFitness, bestPosition] = swoop(pop, nPop, bestFitness, bestPosition, lb, ub, dim, objFun);

        convergenceCurve(t) = bestFitness;
    end
end

% --- Helper Functions --- %

function [pop, bestF, bestX] = select_space(pop, nPop, bestF, bestX, lb, ub, dim, fobj)
    Mean = mean(pop.pos, 1);
    lm = 2;
    for i = 1:nPop
        newPos = bestX + lm*rand(1,dim).*(Mean - pop.pos(i,:));
        newPos = max(min(newPos, ub), lb);
        newCost = fobj(newPos);
        if newCost < pop.cost(i)
            pop.pos(i,:) = newPos;
            pop.cost(i) = newCost;
            if newCost < bestF
                bestX = newPos;
                bestF = newCost;
            end
        end
    end
end

function [pop, bestF, bestX] = search_space(pop, nPop, bestF, bestX, lb, ub, dim, fobj)
    Mean = mean(pop.pos, 1);
    a = 10; R = 1.5;
    A = randperm(nPop);
    pop.pos = pop.pos(A,:);
    pop.cost = pop.cost(A);
    [x, y] = polr(a,R,nPop);
    for i = 1:nPop-1
        Step = pop.pos(i,:) - pop.pos(i+1,:);
        Step1 = pop.pos(i,:) - Mean;
        newPos = pop.pos(i,:) + y(i)*Step + x(i)*Step1;
        newPos = max(min(newPos, ub), lb);
        newCost = fobj(newPos);
        if newCost < pop.cost(i)
            pop.pos(i,:) = newPos;
            pop.cost(i) = newCost;
            if newCost < bestF
                bestX = newPos;
                bestF = newCost;
            end
        end
    end
end

function [pop, bestF, bestX] = swoop(pop, nPop, bestF, bestX, lb, ub, dim, fobj)
    Mean = mean(pop.pos,1);
    a = 10; R = 1.5;
    A = randperm(nPop);
    pop.pos = pop.pos(A,:);
    pop.cost = pop.cost(A);
    [x, y] = swoo_p(a,R,nPop);
    for i = 1:nPop
        Step = pop.pos(i,:) - 2*Mean;
        Step1 = pop.pos(i,:) - 2*bestX;
        newPos = rand(1,dim).*bestX + x(i)*Step + y(i)*Step1;
        newPos = max(min(newPos, ub), lb);
        newCost = fobj(newPos);
        if newCost < pop.cost(i)
            pop.pos(i,:) = newPos;
            pop.cost(i) = newCost;
            if newCost < bestF
                bestX = newPos;
                bestF = newCost;
            end
        end
    end
end

function [xR, yR] = swoo_p(a, R, N)
    th = a*pi*exp(rand(N,1));
    r = th;
    xR = r.*sinh(th);
    yR = r.*cosh(th);
    xR = xR/max(abs(xR));
    yR = yR/max(abs(yR));
end

function [xR, yR] = polr(a, R, N)
    th = a*pi*rand(N,1);
    r = th + R*rand(N,1);
    xR = r.*sin(th);
    yR = r.*cos(th);
    xR = xR/max(abs(xR));
    yR = yR/max(abs(yR));
end


