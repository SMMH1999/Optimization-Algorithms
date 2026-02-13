function [bestFitness, bestPosition, convergenceCurve] = RFO(lb, ub, dim, nPop, maxItr, objFun)
    %% Rupple Fox Optimizer (RFO)
    % Author: Laith Abualigah (adapted)
    %
    % Description:
    %   RFO is a nature-inspired metaheuristic optimizer based on the
    %   behaviors of Rupple Foxes (eyesight, hearing, smell, and social interactions)
    %   for continuous, unconstrained, single-objective optimization problems.
    %
    % Inputs:
    %   lb       - lower bound (scalar or 1 x dim vector)
    %   ub       - upper bound (scalar or 1 x dim vector)
    %   dim      - number of decision variables
    %   nPop     - population size (number of Rupple Foxes)
    %   maxItr   - maximum number of iterations
    %   objFun   - handle to the objective function to minimize
    %
    % Outputs:
    %   bestFitness      - best objective value found
    %   bestPosition     - position of the global best solution
    %   convergenceCurve - vector containing best fitness at each iteration

    %% Initialization
    convergenceCurve = zeros(1, maxItr);
    beta = 1e-10;        % small constant for random walk step control
    vec_flag = [1, -1];  % direction flags for random walk

    % Initialize population
    Positions = initializePopulation(nPop, dim, ub, lb);

    % Evaluate initial fitness
    fitness = zeros(nPop,1);
    for i = 1:nPop
        fitness(i) = objFun(Positions(i,:));
    end

    [bestFitness, idx] = min(fitness);
    pbest = Positions;        % personal best positions
    bestPosition = pbest(idx,:);
    x = pbest;                % current positions

    L = 100;                  % parameter controlling sight/hearing sigmoid

    %% Main Loop
    for ite = 1:maxItr
        % Eyesight (h) and hearing (s) dynamics
        h = 1/(1 + exp((ite - maxItr/2)/L));
        s = 1/(1 + exp((maxItr/2 - ite)/L));
        smell = 0.1 / abs(acos(2 / (1 + exp((maxItr/2 - ite)/100))));

        %% Behavioral updates
        for i = 1:nPop
            % Choose a random fox
            dr = randi(nPop);
            ds = randi(nPop);
            Flag = vec_flag(randi(2));

            if rand >= 0.5  % Day
                if s >= h
                    if rand >= 0.25
                        Xrand = pbest(ds,:);
                        randIdx = rand * randi([1,4]);
                        Positions(i,:) = x(i,:) + randIdx*(Xrand - x(i,:)) + randIdx*(bestPosition - x(i,:));
                        theta = rand*260*pi/360;
                        Positions(i,:) = rotateVector(Positions(i,:), Xrand, theta, dim);
                    else
                        Xrand = pbest(ds,:);
                        Positions(i,:) = Xrand + beta*randn(1,dim)*Flag;
                    end
                else
                    if rand >= 0.75
                        Xrand = pbest(ds,:);
                        randIdx = rand*randi([1,4]);
                        Positions(i,:) = x(i,:) + randIdx*(Xrand - x(i,:)) + randIdx*(bestPosition - x(i,:));
                        theta = rand*150*pi/180;
                        Positions(i,:) = rotateVector(Positions(i,:), Xrand, theta, dim);
                    else
                        Xrand = pbest(ds,:);
                        Positions(i,:) = Xrand + beta*randn(1,dim)*Flag;
                    end
                end
            else % Night
                if h >= s
                    if rand >= 0.25
                        Xrand = pbest(dr,:);
                        randIdx = rand*randi([1,4]);
                        Positions(i,:) = x(i,:) + randIdx*(Xrand - x(i,:)) + randIdx*(bestPosition - x(i,:));
                        theta = rand*150*pi/360;
                        Positions(i,:) = rotateVector(Positions(i,:), Xrand, theta, dim);
                    else
                        Positions(i,:) = pbest(dr,:) + beta*randn(1,dim)*Flag;
                    end
                else
                    if rand >= 0.75
                        Xrand = pbest(dr,:);
                        randIdx = rand*randi([1,4]);
                        Positions(i,:) = x(i,:) + randIdx*(Xrand - x(i,:)) + randIdx*(bestPosition - x(i,:));
                        theta = rand*260*pi/180;
                        Positions(i,:) = rotateVector(Positions(i,:), Xrand, theta, dim);
                    else
                        Positions(i,:) = pbest(dr,:) + beta*randn(1,dim)*Flag;
                    end
                end
            end
        end

        %% Smell behavior
        for i = 1:nPop
            ss = randi(nPop);
            Flag = vec_flag(randi(2));
            if rand >= smell
                Xrand = pbest(ss,:);
                xr = 2 + rand*(4-2);
                eps = abs(4*rand - (rand+rand))/xr;
                Positions(i,:) = x(i,:) + (Xrand - x(i,:))*eps + eps*(bestPosition - x(i,:));
            else
                Positions(i,:) = pbest(ss,:) + beta*randn(1,dim)*Flag;
            end
        end

        %% Animal social behavior
        for i = 1:nPop
            ab = randi(nPop);
            Dista = 2*(bestPosition - Positions(i,:))*rand;
            Distb = 3*(pbest(ab,:) - pbest(i,:))*rand;
            Flag = vec_flag(randi(2));
            if rand >= h
                Positions(i,:) = pbest(ab,:) + beta*randn(1,dim);
            else
                if i == 1
                    Positions(i,:) = Positions(i,:) + Dista + Distb*Flag;
                else
                    xNew = Positions(i,:) + Dista + Distb*Flag;
                    Positions(i,:) = (xNew + Positions(i-1,:))/2;
                end
            end
        end

        %% Worst-case exploration
        if rand > 0.5
            [~, worstIdx] = max(fitness);
            Flag = vec_flag(randi(2));
            Positions(worstIdx,:) = Positions(worstIdx,:) + beta*randn(1,dim)*Flag;
        end

        %% Boundary handling
        Positions = max(min(Positions, ub), lb);

        %% Fitness evaluation & update personal/global best
        for i = 1:nPop
            fitVal = objFun(Positions(i,:));
            if fitVal < fitness(i)
                fitness(i) = fitVal;
                pbest(i,:) = Positions(i,:);
            end
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = pbest(i,:);
            end
        end

        %% Store convergence
        convergenceCurve(ite) = bestFitness;
    end

end

%% Helper: Initialize Population
function Positions = initializePopulation(nPop, dim, ub, lb)
    if isscalar(ub) && isscalar(lb)
        Positions = rand(nPop, dim).*(ub-lb) + lb;
    else
        Positions = zeros(nPop, dim);
        for j = 1:dim
            Positions(:,j) = rand(nPop,1)*(ub(j)-lb(j)) + lb(j);
        end
    end
end

%% Helper: Rotation operation
function v_rot = rotateVector(x, y, theta, dim)
    cent = ceil(dim/2);
    v = [x; y];
    x_center = x(cent);
    y_center = y(cent);
    center = repmat([x_center; y_center], 1, dim);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    vo = R*(v - center) + center;
    v_rot = vo(1,:);
end
