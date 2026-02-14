function [bestFitness, bestPosition, convergenceCurve] = MGO(lb, ub, dim, nPop, maxItr, objFun)
    %% Moss Growth Optimization (MGO)
    % Author: Boli Zheng et al., 2024
    % Journal: Journal of Computational Design and Engineering
    %
    % Description:
    %   Implementation of the Moss Growth Optimization (MGO) algorithm for
    %   minimization problems.
    %
    % Inputs:
    %   lb       - lower bound (scalar or 1xdim vector)
    %   ub       - upper bound (scalar or 1xdim vector)
    %   dim      - number of decision variables
    %   nPop     - population size (number of moss agents)
    %   maxItr   - maximum number of iterations
    %   objFun   - handle of the objective function
    %
    % Outputs:
    %   bestFitness       - best objective function value found
    %   bestPosition      - decision variables corresponding to bestFitness
    %   convergenceCurve  - vector of bestFitness values over iterations

    %% Initialization
    bestFitness = inf; % For minimization
    bestPosition = zeros(1, dim);

    M = initialization(nPop, dim, ub, lb); % Initialize population
    costs = zeros(1, nPop);

    for i = 1:nPop
        costs(i) = objFun(M(i,:));
        if costs(i) < bestFitness
            bestPosition = M(i,:);
            bestFitness = costs(i);
        end
    end

    convergenceCurve = zeros(1, maxItr);

    % Algorithm parameters
    w = 2;
    rec_num = 10;
    divide_num = dim / 4;
    d1 = 0.2;

    newM = zeros(nPop, dim);
    newM_cost = zeros(1, nPop);
    rM = zeros(nPop, dim, rec_num);
    rM_cos = zeros(1, nPop, rec_num);
    rec = 1;

    %% Main Loop
    for it = 1:maxItr
        calPositions = M;
        div_idx = randperm(dim);

        % Divide population and select region with majority
        for j = 1:max(divide_num,1)
            th = bestPosition(div_idx(j));
            index = calPositions(:, div_idx(j)) > th;
            if sum(index) < size(calPositions,1)/2
                index = ~index;
            end
            calPositions = calPositions(index,:);
        end

        D = bestPosition - calPositions;
        D_wind = sum(D,1) / size(calPositions,1);

        beta = size(calPositions,1) / nPop;
        gamma = 1 / sqrt(1 - beta^2);

        step  = w * (rand(size(D_wind))-0.5);
        step2 = 0.1 * w * (rand(size(D_wind))-0.5) * (1 + 0.5*(1 + tanh(beta/gamma)));
        step3 = 0.1 * (rand()-0.5);
        act = actCal(1 ./ (1 + (0.5 - 10*(rand(size(D_wind))))));

        if rec == 1
            rM(:,:,rec) = M;
            rM_cos(1,:,rec) = costs;
            rec = rec + 1;
        end

        % Update each moss agent
        for i = 1:nPop
            newM(i,:) = M(i,:);

            % Spore dispersal search
            if rand() > d1
                newM(i,:) = newM(i,:) + step .* D_wind;
            else
                newM(i,:) = newM(i,:) + step2 .* D_wind;
            end

            % Dual propagation search
            if rand() < 0.8
                if rand() > 0.5
                    newM(i, div_idx(1)) = bestPosition(div_idx(1)) + step3 * D_wind(div_idx(1));
                else
                    newM(i,:) = (1 - act) .* newM(i,:) + act .* bestPosition;
                end
            end

            % Boundary check
            newM(i,:) = max(min(newM(i,:), ub), lb);

            % Evaluate
            newM_cost(i) = objFun(newM(i,:));

            % Record history for cryptobiosis
            rM(i,:,rec) = newM(i,:);
            rM_cos(1,i,rec) = newM_cost(i);

            % Update global best
            if newM_cost(i) < bestFitness
                bestPosition = newM(i,:);
                bestFitness = newM_cost(i);
            end
        end

        rec = rec + 1;

        % Cryptobiosis mechanism
        if rec > rec_num
            [lcost, Iindex] = min(rM_cos, [], 3);
            for i = 1:nPop
                M(i,:) = rM(i,:,Iindex(i));
            end
            costs = lcost;
            rec = 1;
        else
            M = newM;
            costs = newM_cost;
        end

        convergenceCurve(it) = bestFitness;
    end

end

%% Helper function: action calculation
function [act] = actCal(X)
    act = X;
    act(act >= 0.5) = 1;
    act(act < 0.5) = 0;
end

%% Helper function: initialize population
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Boundary_no = numel(ub);
    Positions = zeros(SearchAgents_no, dim);

    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) * (ub - lb) + lb;
    else
        for i = 1:dim
            Positions(:,i) = rand(SearchAgents_no,1) * (ub(i)-lb(i)) + lb(i);
        end
    end
end
