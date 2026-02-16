function [bestFitness, bestPosition, convergenceCurve] = FHO(lb, ub, dim, nPop, maxItr, objFun)
    %__________________________________________________________________________
    % Fire Hawk Optimizer (FHO) - Version 1.0
    %
    % Author: Mahdi Azizi, Siamak TalatAhari, Amir H. Gandomi
    % Main Source: Fire Hawk Optimizer: a novel metaheuristic algorithm
    %              https://link.springer.com/article/10.1007/s10462-022-10173-w
    %
    % Description:
    %   Metaheuristic optimization algorithm inspired by Fire Hawks hunting behavior.
    %
    % Inputs:
    %   lb        - Lower bound (scalar or 1xDim vector)
    %   ub        - Upper bound (scalar or 1xDim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size
    %   maxItr    - Maximum number of iterations
    %   objFun    - Handle to objective function: fitness = objFun(position)
    %
    % Outputs:
    %   bestFitness       - Best fitness value found
    %   bestPosition      - Decision vector corresponding to bestFitness
    %   convergenceCurve  - Vector of bestFitness at each iteration
    %__________________________________________________________________________

    %% Ensure bounds are row vectors
    if isscalar(lb), lb = lb*ones(1,dim); end
    if isscalar(ub), ub = ub*ones(1,dim); end

    %% Initialization
    Pop = unifrnd(lb, ub, nPop, dim);
    Cost = zeros(nPop,1);

    for i = 1:nPop
        Cost(i) = objFun(Pop(i,:));
    end

    [Cost, idx] = sort(Cost);
    Pop = Pop(idx,:);
    BestPop = Pop(1,:);
    SP = mean(Pop,1);

    HN = randi([1 ceil(nPop/5)],1,1);      % Number of Fire Hawks
    FHPops = Pop(1:HN,:);
    Pop2 = Pop(HN+1:end,:);

    PopNew = cell(HN,1);
    Pop2_temp = Pop2;
    for i = 1:HN
        nPop2 = size(Pop2_temp,1);
        if nPop2 < HN, break; end
        Dist = vecnorm(FHPops(i,:) - Pop2_temp, 2, 2);
        [~, b] = sort(Dist);
        alfa = randi(nPop2);
        PopNew{i} = Pop2_temp(b(1:alfa),:);
        Pop2_temp(b(1:alfa),:) = [];
        if isempty(Pop2_temp), break; end
    end
    if ~isempty(Pop2_temp)
        PopNew{end} = [PopNew{end}; Pop2_temp];
    end

    GB = Cost(1);           % Global best fitness
    BestPos = BestPop;      % Global best position
    convergenceCurve = [];

    %% Main Loop
    for Iter = 1:maxItr
        PopTot = [];

        for i = 1:length(PopNew)
            PR = PopNew{i};
            FHl = FHPops(i,:);
            SPl = mean(PR,1);

            % Update Fire Hawk
            Ir = rand(1,2);
            FHnear = FHPops(randi(HN),:);
            FHl_new = FHl + (Ir(1)*GB - Ir(2)*FHnear);
            FHl_new = min(max(FHl_new, lb), ub);
            PopTot = [PopTot; FHl_new];

            % Update Prey
            for q = 1:size(PR,1)
                Ir = rand(1,2);
                PRq_new1 = PR(q,:) + (Ir(1)*FHl - Ir(2)*SPl);
                PRq_new1 = min(max(PRq_new1, lb), ub);
                PopTot = [PopTot; PRq_new1];

                Ir = rand(1,2);
                FHAlter = FHPops(randi(HN),:);
                PRq_new2 = PR(q,:) + (Ir(1)*FHAlter - Ir(2)*SP);
                PRq_new2 = min(max(PRq_new2, lb), ub);
                PopTot = [PopTot; PRq_new2];
            end
        end

        % Evaluate new population
        CostTot = zeros(size(PopTot,1),1);
        for i = 1:size(PopTot,1)
            CostTot(i) = objFun(PopTot(i,:));
        end

        % Sort and select top nPop
        [CostTot, idx] = sort(CostTot);
        PopTot = PopTot(idx,:);
        Pop = PopTot(1:nPop,:);

        % Update Fire Hawks and Prey
        HN = randi([1 ceil(nPop/5)],1,1);
        BestPop = Pop(1,:);
        SP = mean(Pop,1);
        FHPops = Pop(1:HN,:);
        Pop2 = Pop(HN+1:end,:);

        PopNew = cell(HN,1);
        Pop2_temp = Pop2;
        for i = 1:HN
            nPop2 = size(Pop2_temp,1);
            if nPop2 < HN, break; end
            Dist = vecnorm(FHPops(i,:) - Pop2_temp, 2, 2);
            [~, b] = sort(Dist);
            alfa = randi(nPop2);
            PopNew{i} = Pop2_temp(b(1:alfa),:);
            Pop2_temp(b(1:alfa),:) = [];
            if isempty(Pop2_temp), break; end
        end
        if ~isempty(Pop2_temp)
            PopNew{end} = [PopNew{end}; Pop2_temp];
        end

        % Update global best
        if CostTot(1) < GB
            BestPos = Pop(1,:);
        end
        GB = min(GB, CostTot(1));

        convergenceCurve(Iter,1) = GB;
    end

    bestFitness = GB;
    bestPosition = BestPos;

end
