function [bestFitness, bestPosition, convergenceCurve] = PO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Political Optimizer (PO)
    % =========================================================================
    % Author: Qamar Askari
    % Source Paper:
    % Askari, Q., Younas, I., & Saeed, M. (2020).
    % Political Optimizer: A novel socio-inspired meta-heuristic for global optimization.
    % Knowledge-Based Systems.
    %
    % -------------------------------------------------------------------------
    % Description:
    % Political Optimizer (PO) is a socio-inspired metaheuristic algorithm
    % based on political processes including election, party switching,
    % government formation, and parliamentarism.
    %
    % -------------------------------------------------------------------------
    % Inputs:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Number of decision variables
    % nPop      : Population size (must be divisible into parties × areas)
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % -------------------------------------------------------------------------
    % Outputs:
    % bestFitness      : Best objective function value found
    % bestPosition     : Best decision vector found
    % convergenceCurve : Best fitness value at each iteration
    %
    % -------------------------------------------------------------------------
    % Tunable Parameters (Internal):
    % parties : Number of political parties (default = 5)
    % areas   : Members per party (nPop/parties)
    % lambda  : Party switching rate control parameter (default = 1)
    %
    % =========================================================================

    %% Parameter setting
    parties = 5;
    areas   = floor(nPop/parties);
    nPop    = parties * areas;
    lambda  = 1;

    if numel(lb) == 1
        lb = lb * ones(1, dim);
    end
    if numel(ub) == 1
        ub = ub * ones(1, dim);
    end

    %% Initialization
    Positions = rand(nPop, dim) .* (ub - lb) + lb;
    prevPositions = Positions;
    auxPositions  = Positions;

    fitness       = zeros(nPop,1);
    prevFitness   = inf(nPop,1);
    auxFitness    = inf(nPop,1);

    bestFitness   = inf;
    bestPosition  = zeros(1,dim);
    convergenceCurve = zeros(1,maxItr);

    %% Initial Election (Fitness Evaluation)
    for i = 1:nPop
        fitness(i) = objFun(Positions(i,:));
        if fitness(i) < bestFitness
            bestFitness  = fitness(i);
            bestPosition = Positions(i,:);
        end
    end

    %% Government Formation
    [aWinnerInd, aWinners, pLeaderInd, pLeaders] = governmentFormation();

    %% Main Loop
    for t = 1:maxItr

        prevFitness   = auxFitness;
        prevPositions = auxPositions;
        auxFitness    = fitness;
        auxPositions  = Positions;

        %% Election Campaign
        for whichMethod = 1:2
            for a = 1:areas
                for p = 1:parties
                    i = (p-1)*areas + a;
                    for j = 1:dim

                        if whichMethod == 1
                            center = pLeaders(p,j);
                        else
                            center = aWinners(a,j);
                        end

                        if prevFitness(i) >= fitness(i)
                            if (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) <= center) || ...
                                    (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) >= center)
                                radius = center - Positions(i,j);
                                Positions(i,j) = center + rand()*radius;
                            else
                                radius = abs(Positions(i,j) - center);
                                Positions(i,j) = center + (2*rand()-1)*radius;
                            end
                        else
                            if (prevPositions(i,j) <= Positions(i,j) && Positions(i,j) <= center) || ...
                                    (prevPositions(i,j) >= Positions(i,j) && Positions(i,j) >= center)
                                radius = abs(Positions(i,j) - center);
                                Positions(i,j) = center + (2*rand()-1)*radius;
                            else
                                radius = Positions(i,j) - prevPositions(i,j);
                                Positions(i,j) = prevPositions(i,j) + rand()*radius;
                            end
                        end

                    end
                end
            end
        end

        %% Party Switching
        psr = (1 - t/maxItr) * lambda;
        for p = 1:parties
            for a = 1:areas
                fromInd = (p-1)*areas + a;
                if rand < psr
                    toParty = randi(parties);
                    while toParty == p
                        toParty = randi(parties);
                    end
                    toStart = (toParty-1)*areas + 1;
                    toEnd   = toStart + areas - 1;
                    [~, worstIdx] = max(fitness(toStart:toEnd));
                    toInd = toStart + worstIdx - 1;

                    temp = Positions(toInd,:);
                    Positions(toInd,:) = Positions(fromInd,:);
                    Positions(fromInd,:) = temp;

                    temp = fitness(toInd);
                    fitness(toInd) = fitness(fromInd);
                    fitness(fromInd) = temp;
                end
            end
        end

        %% Election Phase (Evaluation + Bound Control)
        for i = 1:nPop
            Positions(i,:) = max(Positions(i,:), lb);
            Positions(i,:) = min(Positions(i,:), ub);

            fitness(i) = objFun(Positions(i,:));
            if fitness(i) < bestFitness
                bestFitness  = fitness(i);
                bestPosition = Positions(i,:);
            end
        end

        %% Government Formation
        [aWinnerInd, aWinners, pLeaderInd, pLeaders] = governmentFormation();

        %% Parliamentarism
        for a = 1:areas
            i = aWinnerInd(a);
            newAWinner = aWinners(a,:);

            toa = randi(areas);
            while toa == a
                toa = randi(areas);
            end

            toAWinner = aWinners(toa,:);
            for j = 1:dim
                distance = abs(toAWinner(j) - newAWinner(j));
                newAWinner(j) = toAWinner(j) + (2*rand()-1)*distance;
            end

            newFit = objFun(newAWinner);
            if newFit < fitness(i)
                Positions(i,:) = newAWinner;
                fitness(i) = newFit;
                aWinners(a,:) = newAWinner;

                if newFit < bestFitness
                    bestFitness  = newFit;
                    bestPosition = newAWinner;
                end
            end
        end

        convergenceCurve(t) = bestFitness;
    end

    %% ================= Local Function =================
    function [aWinnerInd, aWinners, pLeaderInd, pLeaders] = governmentFormation()

        aWinnerInd = zeros(areas,1);
        aWinners   = zeros(areas,dim);

        for a = 1:areas
            [~, idx] = min(fitness(a:areas:nPop));
            aWinnerInd(a) = (idx-1)*areas + a;
            aWinners(a,:) = Positions(aWinnerInd(a),:);
        end

        pLeaderInd = zeros(parties,1);
        pLeaders   = zeros(parties,dim);

        for p = 1:parties
            st = (p-1)*areas + 1;
            en = st + areas - 1;
            [~, idx] = min(fitness(st:en));
            pLeaderInd(p) = st + idx - 1;
            pLeaders(p,:) = Positions(pLeaderInd(p),:);
        end

    end

end
