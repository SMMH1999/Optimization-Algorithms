function [bestFitness, bestPosition, convergenceCurve] = SFOA(lb, ub, dim, nPop, maxItr, objFun)
    % Starfish Optimization Algorithm (SFOA)
    % Abbreviation: SFOA
    % Author: Dr. Changting Zhong
    % Paper: Zhong et al., Neural Computing and Applications, 2025, 37: 3641-3683
    %
    % Description:
    %   SFOA is a bio-inspired metaheuristic algorithm based on starfish behaviors
    %   for global optimization. It balances exploration and exploitation
    %   through probabilistic updates and multiple “arms” interactions.
    %
    % Inputs:
    %   lb       : 1xdim vector or scalar, lower bounds of variables
    %   ub       : 1xdim vector or scalar, upper bounds of variables
    %   dim      : integer, number of decision variables
    %   nPop     : integer, population size
    %   maxItr   : integer, maximum number of iterations
    %   objFun   : function handle, objective function to minimize
    %
    % Outputs:
    %   bestFitness      : scalar, best objective value found
    %   bestPosition     : 1xdim vector, position corresponding to bestFitness
    %   convergenceCurve : 1xmaxItr vector, bestFitness at each iteration
    %
    % Tunable Parameters:
    %   GP : exploration probability (default 0.5)

    %% Parameter
    GP = 0.5;  % exploration probability

    %% Boundary handling
    if isscalar(ub)
        lb = lb*ones(1,dim);
        ub = ub*ones(1,dim);
    end

    %% Initialization
    X = rand(nPop,dim).*(ub-lb) + lb;      % population initialization
    Fitness = zeros(1,nPop);

    for i = 1:nPop
        Fitness(i) = objFun(X(i,:));
    end

    [bestFitness, idx] = min(Fitness);      % global best fitness
    bestPosition = X(idx,:);                % global best position

    newX = zeros(nPop,dim);
    convergenceCurve = zeros(1,maxItr);

    %% Main loop
    for T = 1:maxItr
        theta = pi/2 * T / maxItr;
        tEO = (maxItr-T)/maxItr * cos(theta);

        if rand < GP   % Exploration
            for i = 1:nPop
                if dim > 5
                    jp = randperm(dim,5);
                    for j = 1:5
                        pm = (2*rand-1)*pi;
                        if rand < GP
                            newX(i,jp(j)) = X(i,jp(j)) + pm*(bestPosition(jp(j))-X(i,jp(j)))*cos(theta);
                        else
                            newX(i,jp(j)) = X(i,jp(j)) - pm*(bestPosition(jp(j))-X(i,jp(j)))*sin(theta);
                        end
                        % Boundary check
                        if newX(i,jp(j)) > ub(jp(j)) || newX(i,jp(j)) < lb(jp(j))
                            newX(i,jp(j)) = X(i,jp(j));
                        end
                    end
                else
                    jp = ceil(dim*rand);
                    im = randperm(nPop,2);
                    rand1 = 2*rand-1; rand2 = 2*rand-1;
                    newX(i,jp) = tEO*X(i,jp) + rand1*(X(im(1),jp)-X(i,jp)) + rand2*(X(im(2),jp)-X(i,jp));
                    if newX(i,jp) > ub(jp) || newX(i,jp) < lb(jp)
                        newX(i,jp) = X(i,jp);
                    end
                end
                newX(i,:) = max(min(newX(i,:),ub),lb);
            end
        else  % Exploitation
            df = randperm(nPop,5);
            dm = bestPosition - X(df,:);
            for i = 1:nPop
                r1 = rand; r2 = rand;
                kp = randperm(5,2);
                newX(i,:) = X(i,:) + r1*dm(kp(1),:) + r2*dm(kp(2),:);
                if i == nPop
                    newX(i,:) = exp(-T*nPop/maxItr) * X(i,:); % regeneration
                end
                newX(i,:) = max(min(newX(i,:),ub),lb);
            end
        end

        % Fitness evaluation and update
        for i = 1:nPop
            fit = objFun(newX(i,:));
            if fit < Fitness(i)
                Fitness(i) = fit;
                X(i,:) = newX(i,:);
                if fit < bestFitness
                    bestFitness = fit;
                    bestPosition = X(i,:);
                end
            end
        end

        convergenceCurve(T) = bestFitness;
    end

end
