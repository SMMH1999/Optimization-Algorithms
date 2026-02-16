function [bestFitness, bestPosition, convergenceCurve] = FOA(lb, ub, dim, nPop, maxItr, objFun)
    % FOA - Forest Optimization Algorithm (safe, no resizing in loops)
    %
    % Inputs:
    %   lb, ub       : scalar or vector bounds (1 x dim)
    %   dim          : number of variables
    %   nPop         : number of trees
    %   maxItr       : maximum iterations
    %   objFun       : function handle for fitness evaluation
    %
    % Outputs:
    %   bestFitness      : best fitness found
    %   bestPosition     : position of best fitness
    %   convergenceCurve : best fitness per iteration

    lb = lb(:)'; ub = ub(:)';  % ensure row vectors

    % Algorithm parameters
    dx = abs(ub - lb)/5;
    lifeTime = ceil(nPop/2);
    transferRate = 15;
    LSC = max(1, floor(2*dim/10));
    GSC = max(1, floor(1*dim/10));

    % Initialize Population and helper arrays
    Population = rand(nPop, dim) .* (ub-lb) + lb;
    fitness = zeros(nPop,1);
    age = zeros(nPop,1);

    for i = 1:nPop
        fitness(i) = objFun(Population(i,:));
    end

    [bestFitness, idx] = min(fitness);
    bestPosition = Population(idx,:);

    convergenceCurve = zeros(1,maxItr);

    % Preallocate maximum possible size for new trees per iteration
    maxNew = nPop * LSC;  % upper bound for local seeding
    newPop = zeros(maxNew, dim);
    newFitness = zeros(maxNew,1);
    newAge = zeros(maxNew,1);

    for iter = 1:maxItr
        % ---------------- Local Seeding ----------------
        newCount = 0;
        for i = 1:nPop
            if age(i) == 0
                varsToChange = randperm(dim, LSC);
                for v = varsToChange
                    newCount = newCount + 1;
                    tempVar = Population(i,:);
                    tempVar(v) = tempVar(v) + (rand*2-1)*dx(v);
                    tempVar(v) = min(max(tempVar(v), lb(v)), ub(v));
                    newPop(newCount,:) = tempVar;
                    newFitness(newCount) = objFun(tempVar);
                    newAge(newCount) = 0;
                end
            end
        end

        % Only keep the filled part of newPop
        if newCount > 0
            Population = [Population; newPop(1:newCount,:)];
            fitness = [fitness; newFitness(1:newCount)];
            age = [age + 1; newAge(1:newCount)];
        else
            age = age + 1;
        end

        % ---------------- Population Limiting ----------------
        aliveMask = age < lifeTime;
        candidateMask = ~aliveMask;

        candidatePop = Population(candidateMask,:);
        candidateFit = fitness(candidateMask);
        candidateAge = age(candidateMask);

        Population = Population(aliveMask,:);
        fitness = fitness(aliveMask);
        age = age(aliveMask);

        % Keep top nPop trees
        if size(Population,1) > nPop
            [~, idxSort] = sort(fitness);
            extraIdx = idxSort(nPop+1:end);

            candidatePop = [candidatePop; Population(extraIdx,:)];
            candidateFit = [candidateFit; fitness(extraIdx)];
            candidateAge = [candidateAge; age(extraIdx)];

            Population(extraIdx,:) = [];
            fitness(extraIdx) = [];
            age(extraIdx) = [];
        end

        % ---------------- Global Seeding ----------------
        if ~isempty(candidatePop)
            RevoSize = max(1,floor(size(candidatePop,1)*transferRate/100));
            idxVars = randi(dim, 1, GSC);
            for w = 1:RevoSize
                tempVar = candidatePop(w,:);
                tempVar(idxVars) = rand(1,GSC).*(ub(idxVars)-lb(idxVars)) + lb(idxVars);
                candidatePop(w,:) = tempVar;
                candidateFit(w) = objFun(tempVar);
            end

            % Append global seeds safely
            Population = [Population; candidatePop(1:RevoSize,:)];
            fitness = [fitness; candidateFit(1:RevoSize)];
            age = [age; zeros(RevoSize,1)];
        end

        % ---------------- Update Best ----------------
        [minFit, minIdx] = min(fitness);
        if minFit < bestFitness
            bestFitness = minFit;
            bestPosition = Population(minIdx,:);
        end

        convergenceCurve(iter) = bestFitness;
    end

end
