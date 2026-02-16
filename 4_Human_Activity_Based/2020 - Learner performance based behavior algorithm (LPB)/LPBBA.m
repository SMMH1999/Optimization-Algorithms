function [bestFitness, bestPosition, convergenceCurve] = LPBBA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Learner Performance Based Behavior Algorithm (LPBBA)
    % -------------------------------------------------------------------------
    % Author: Chnoor M. Rahman
    % Citation:
    % Rahman, C. and Rashid, T., 2020. A new evolutionary algorithm:
    % Learner performance based behavior algorithm.
    % Egyptian Informatics Journal.
    % DOI: https://doi.org/10.1016/j.eij.2020.08.003
    %
    % Description:
    % LPBBA is an evolutionary optimization algorithm inspired by learner
    % performance behavior. The population is divided into different
    % performance-based groups and reproduction is performed accordingly.
    %
    % Inputs:
    %   lb        - Lower bound (scalar or 1×dim vector)
    %   ub        - Upper bound (scalar or 1×dim vector)
    %   dim       - Number of decision variables
    %   nPop      - Population size
    %   maxItr    - Maximum number of iterations
    %   objFun    - Objective function handle (minimization)
    %
    % Outputs:
    %   bestFitness      - Best fitness value found
    %   bestPosition     - Best solution vector
    %   convergenceCurve - Best fitness at each iteration
    %
    % Tunable Parameters:
    %   pc    - Crossover percentage
    %   pm    - Mutation percentage
    %   gamma - Crossover expansion factor
    %   mu    - Mutation rate
    %   dp    - Division probability
    %   beta  - Selection pressure coefficient
    % =========================================================================

    %% Parameter Settings
    pc    = 0.6;
    pm    = 0.3;
    gamma = 0.05;
    mu    = 0.03;
    dp    = 0.5;
    beta  = 8;

    nc = 2 * round(pc * nPop / 2);
    nm = round(pm * nPop);

    %% Bound Handling
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    %% Initialization
    empty.Position = [];
    empty.Cost     = [];

    pop = repmat(empty, nPop, 1);

    for i = 1:nPop
        pop(i).Position = lb + rand(1, dim) .* (ub - lb);
        pop(i).Cost = objFun(pop(i).Position);
    end

    [~, idx] = sort([pop.Cost]);
    pop = pop(idx);

    bestSol = pop(1);
    worstCost = pop(end).Cost;

    convergenceCurve = zeros(maxItr,1);

    %% Main Loop
    for t = 1:maxItr

        Costs = [pop.Cost];
        P = exp(-beta * Costs / worstCost);
        P = P / sum(P);

        [fPop, wPop, gPop, pPop] = dividePopulation(pop, dp);

        val = floor(findValue(numel(fPop)) / 2);

        popc2 = repmat(empty, ceil(0.25*nPop), 2);
        popc2Index = 1;
        ncCheck = 1;

        fIndexUsed = [];

        % First Reproduction Phase
        for c = 1:ceil(val/2)
            i1 = parentSelection(fPop, fIndexUsed);
            fIndexUsed(end+1) = i1;
            i2 = parentSelection(fPop, fIndexUsed);
            fIndexUsed(end+1) = i2;

            parent1 = fPop(i1);
            parent2 = fPop(i2);

            if ncCheck < nc
                popc2(popc2Index,1) = parent1;
                popc2(popc2Index,2) = parent2;

                [child1, child2] = Crossover(parent1.Position,parent2.Position,gamma,lb,ub);

                popc2(popc2Index+1,1).Position = child1;
                popc2(popc2Index+1,2).Position = child2;
                popc2(popc2Index+1,1).Cost = objFun(child1);
                popc2(popc2Index+1,2).Cost = objFun(child2);

                ncCheck = ncCheck + 1;
                popc2Index = popc2Index + 2;
            else
                popc2(popc2Index,1) = parent1;
                popc2(popc2Index,2) = parent2;
                popc2Index = popc2Index + 1;
            end
        end

        popc2 = popc2(:);

        % Mutation
        for k = 1:nm
            i = randi(numel(popc2));
            popc2(i).Position = Mutate(popc2(i).Position,mu,lb,ub);
            popc2(i).Cost = objFun(popc2(i).Position);
        end

        % Merge
        pop = [popc2; pop];
        [~, idx] = sort([pop.Cost]);
        pop = pop(idx);

        worstCost = max(worstCost, pop(end).Cost);
        pop = pop(1:nPop);

        bestSol = pop(1);
        convergenceCurve(t) = bestSol.Cost;
    end

    bestFitness  = bestSol.Cost;
    bestPosition = bestSol.Position;
end
%% ====================== Local Functions ================================

function i1 = parentSelection(fPop, usedIdx)
    i1 = randi(numel(fPop));
    while ismember(i1, usedIdx)
        i1 = randi(numel(fPop));
    end
end

function y = Mutate(x, mu, VarMin, VarMax)
    nVar = numel(x);
    nmu = ceil(mu*nVar);
    j = randperm(nVar, nmu);
    sigma = 0.1*(VarMax-VarMin);
    y = x;
    y(j) = x(j) + sigma(j).*randn(size(j));
    y = max(y, VarMin);
    y = min(y, VarMax);
end

function [y1, y2] = Crossover(x1,x2,gamma,VarMin,VarMax)
    alpha = -gamma + (1+2*gamma).*rand(size(x1));
    y1 = alpha.*x1 + (1-alpha).*x2;
    y2 = alpha.*x2 + (1-alpha).*x1;
    y1 = max(min(y1,VarMax),VarMin);
    y2 = max(min(y2,VarMax),VarMin);
end

function [fPop, wPop, gPop, pPop] = dividePopulation(pop, divProbability)
    nPop = numel(pop);
    sr = floor(nPop * divProbability);
    r = randperm(nPop, sr);
    rPop = pop(r);

    [~, idx] = sort([rPop.Cost]);
    sortedrPop = rPop(idx);

    size1 = floor(sr * divProbability);
    fPop = sortedrPop(1:size1);

    worstOf1 = sortedrPop(end);
    bestIndiv = sortedrPop(1);

    betterPop = pop([pop.Cost] < worstOf1.Cost);
    wPop      = pop([pop.Cost] >= worstOf1.Cost);

    pPop = betterPop([betterPop.Cost] < bestIndiv.Cost);
    gPop = betterPop([betterPop.Cost] >= bestIndiv.Cost);
end

function result = findValue(val)
    if mod(val,4)==0
        result = val;
    else
        result = val + 2;
    end
end


