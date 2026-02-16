function [bestFitness, bestPosition, convergenceCurve] = SPO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Stochastic Paint Optimizer (SPO)
    % =========================================================================
    % Authors: Prof. Ali Kaveh, Nima Khodadadi, Simak Talatahari
    % Reference:
    % Kaveh, A., Talatahari, S., & Khodadadi, N. (2020).
    % Stochastic Paint Optimizer: Theory and Application in Civil Engineering.
    % Engineering with Computers.
    %
    % -------------------------------------------------------------------------
    % Description:
    % SPO is a population-based stochastic optimization algorithm inspired by
    % paint combination strategies. The population (colors) is divided into
    % three ranked groups, and new candidate solutions are generated using:
    %   1) Complement Combination
    %   2) Analog Combination
    %   3) Triangle Combination
    %   4) Rectangle Combination
    % The best individuals survive to the next iteration.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    %   lb        : Lower bound (scalar or 1×dim vector)
    %   ub        : Upper bound (scalar or 1×dim vector)
    %   dim       : Number of decision variables (dimension)
    %   nPop      : Population size (number of colors)
    %   maxItr    : Maximum number of iterations
    %   objFun    : Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      : Best fitness value found
    %   bestPosition     : Best solution vector found
    %   convergenceCurve : Best fitness at each iteration (1×maxItr)
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters:
    %   nPop   : Population size
    %   maxItr : Number of iterations
    %
    % =========================================================================

    %% ------------------------- Initialization -------------------------------

    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end

    % Group sizes
    N1 = floor(nPop/3);
    N2 = floor(nPop/3);
    N3 = nPop - N1 - N2;

    % Initialize population
    Colors = rand(nPop, dim) .* (ub - lb) + lb;
    Fun_eval = zeros(nPop, 1);

    for i = 1:nPop
        Fun_eval(i) = objFun(Colors(i,:));
    end

    % Initialize global best
    [bestFitness, idxBest] = min(Fun_eval);
    bestPosition = Colors(idxBest,:);
    convergenceCurve = zeros(1, maxItr);

    %% --------------------------- Main Loop ----------------------------------

    for t = 1:maxItr

        for ind = 1:nPop

            % Sort population
            [Fun_eval, sortIdx] = sort(Fun_eval);
            Colors = Colors(sortIdx,:);

            % Divide into three groups
            Group1 = Colors(1:N1,:);
            Group2 = Colors(N1+1:N1+N2,:);
            Group3 = Colors(N1+N2+1:nPop,:);

            NewColors = zeros(4, dim);
            Fun_evalNew = zeros(4,1);

            %% Complement Combination
            id1 = randi(N1);
            id2 = randi(N3);
            NewColors(1,:) = Colors(ind,:) + rand(1,dim) .* ...
                (Group1(id1,:) - Group3(id2,:));

            %% Analog Combination
            if ind <= N1
                id = randi(N1,2,1);
                AnalogGroup = Group1;
            elseif ind <= N1+N2
                id = randi(N2,2,1);
                AnalogGroup = Group2;
            else
                id = randi(N3,2,1);
                AnalogGroup = Group3;
            end
            NewColors(2,:) = Colors(ind,:) + rand(1,dim) .* ...
                (AnalogGroup(id(2),:) - AnalogGroup(id(1),:));

            %% Triangle Combination
            id1 = randi(N1);
            id2 = randi(N2);
            id3 = randi(N3);
            NewColors(3,:) = Colors(ind,:) + rand(1,dim) .* ...
                (Group1(id1,:) + Group2(id2,:) + Group3(id3,:)) / 3;

            %% Rectangle Combination
            id1 = randi(N1);
            id2 = randi(N2);
            id3 = randi(N3);
            id4 = randi(nPop);
            NewColors(4,:) = Colors(ind,:) + ...
                (rand(1,dim).*Group1(id1,:) + ...
                rand(1,dim).*Group2(id2,:) + ...
                rand(1,dim).*Group3(id3,:) + ...
                rand(1,dim).*Colors(id4,:)) / 4;

            %% Boundary Control and Evaluation
            for k = 1:4
                NewColors(k,:) = min(max(NewColors(k,:), lb), ub);
                Fun_evalNew(k) = objFun(NewColors(k,:));
            end

            % Merge populations
            Colors = [Colors; NewColors];
            Fun_eval = [Fun_eval; Fun_evalNew];

        end

        %% Survival Selection
        [Fun_eval, sortIdx] = sort(Fun_eval);
        Colors = Colors(sortIdx,:);

        Colors = Colors(1:nPop,:);
        Fun_eval = Fun_eval(1:nPop);

        %% Update Global Best
        if Fun_eval(1) < bestFitness
            bestFitness = Fun_eval(1);
            bestPosition = Colors(1,:);
        end

        %% Convergence History
        convergenceCurve(t) = bestFitness;

    end

end
