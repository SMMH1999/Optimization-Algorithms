function [bestFitness, bestPosition, convergenceCurve] = DMOA(lb, ub, dim, nPop, maxItr, objFun)
    %_______________________________________________________________________________________%
    %  Dwarf Mongoose Optimization Algorithm (DMOA)                                         %
    %                                                                                       %
    %  Developed originally in MATLAB R2015a (7.13)                                         %
    %  Original Authors: Jeffrey O. Agushaka, Absalom E. Ezugwu, Laith Abualigah            %
    %                                                                                       %
    %  Refactored to Instructor / Benchmark format                                          %
    %                                                                                       %
    %  Algorithm Description:                                                               %
    %  DMOA is a nature-inspired metaheuristic based on the social behavior of dwarf       %
    %  mongooses. The population is divided into Alpha group, Scouts, and Babysitters.     %
    %  The algorithm uses vocalization-driven position updates, sleeping mould behavior,   %
    %  and babysitter replacement to balance exploration and exploitation.                 %
    %                                                                                       %
    %  Inputs:                                                                              %
    %    lb        - Lower bound (scalar or 1×dim vector)                                   %
    %    ub        - Upper bound (scalar or 1×dim vector)                                   %
    %    dim       - Number of decision variables                                           %
    %    nPop      - Population size                                                        %
    %    maxItr    - Maximum number of iterations                                           %
    %    objFun    - Objective function handle (minimization)                               %
    %                                                                                       %
    %  Outputs:                                                                             %
    %    bestFitness      - Best objective value found                                      %
    %    bestPosition     - Best decision vector found (1×dim)                              %
    %    convergenceCurve - Best fitness value at each iteration (maxItr×1)                %
    %                                                                                       %
    %  Internal Parameters:                                                                 %
    %    nBabysitter - Number of babysitters (3)                                            %
    %    peep        - Alpha female vocalization coefficient (2)                            %
    %_______________________________________________________________________________________%

    %% Bound Handling
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    %% Parameters
    VarSize = [1 dim];
    nBabysitter = 3;
    nAlphaGroup = nPop - nBabysitter;
    nScout = nAlphaGroup;
    L = round(0.6 * dim * nBabysitter);
    peep = 2;

    %% Initialization
    empty_mongoose.Position = [];
    empty_mongoose.Cost = [];

    pop = repmat(empty_mongoose, nAlphaGroup, 1);

    BestSol.Cost = inf;
    BestSol.Position = zeros(1, dim);

    tau = inf;
    Iter = 1;
    sm = inf(nAlphaGroup,1);
    C = zeros(nAlphaGroup,1);
    CF = (1 - Iter/maxItr)^(2*Iter/maxItr);

    for i = 1:nAlphaGroup
        pop(i).Position = lb + rand(VarSize).*(ub - lb);
        pop(i).Position = max(pop(i).Position, lb);
        pop(i).Position = min(pop(i).Position, ub);
        pop(i).Cost = objFun(pop(i).Position);

        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end

    convergenceCurve = zeros(maxItr,1);

    %% Main Loop
    for t = 1:maxItr

        %% Alpha Group
        F = zeros(nAlphaGroup,1);
        MeanCost = mean([pop.Cost]);

        for i = 1:nAlphaGroup
            F(i) = exp(-pop(i).Cost / MeanCost);
        end

        P = F / sum(F);

        for m = 1:nAlphaGroup

            i = RouletteWheelSelection(P);
            K = [1:i-1 i+1:nAlphaGroup];
            k = K(randi([1 numel(K)]));

            phi = (peep/2) * (-1 + 2*rand(VarSize));

            newpop.Position = pop(i).Position + phi .* (pop(i).Position - pop(k).Position);
            newpop.Position = max(newpop.Position, lb);
            newpop.Position = min(newpop.Position, ub);

            newpop.Cost = objFun(newpop.Position);

            if newpop.Cost <= pop(i).Cost
                pop(i) = newpop;
            else
                C(i) = C(i) + 1;
            end
        end

        %% Scout Group
        for i = 1:nScout

            K = [1:i-1 i+1:nAlphaGroup];
            k = K(randi([1 numel(K)]));

            phi = (peep/2) * (-1 + 2*rand(VarSize));

            newpop.Position = pop(i).Position + phi .* (pop(i).Position - pop(k).Position);
            newpop.Position = max(newpop.Position, lb);
            newpop.Position = min(newpop.Position, ub);

            newpop.Cost = objFun(newpop.Position);

            sm(i) = (newpop.Cost - pop(i).Cost) / max(newpop.Cost, pop(i).Cost);

            if newpop.Cost <= pop(i).Cost
                pop(i) = newpop;
            else
                C(i) = C(i) + 1;
            end
        end

        %% Babysitters
        for i = 1:nBabysitter
            if C(i) >= L
                pop(i).Position = lb + rand(VarSize).*(ub - lb);
                pop(i).Position = max(pop(i).Position, lb);
                pop(i).Position = min(pop(i).Position, ub);
                pop(i).Cost = objFun(pop(i).Position);
                C(i) = 0;
            end
        end

        %% Update Best
        for i = 1:nAlphaGroup
            if pop(i).Cost <= BestSol.Cost
                BestSol = pop(i);
            end
        end

        %% Next Mongoose Position
        newtau = mean(sm);
        for i = 1:nScout

            M = (pop(i).Position .* sm(i)) ./ pop(i).Position;

            if newtau > tau
                newPosition = pop(i).Position - CF * phi .* rand .* (pop(i).Position - M);
            else
                newPosition = pop(i).Position + CF * phi .* rand .* (pop(i).Position - M);
            end

            newPosition = max(newPosition, lb);
            newPosition = min(newPosition, ub);

            tau = newtau;
        end

        %% Final Best Update
        for i = 1:nAlphaGroup
            if pop(i).Cost <= BestSol.Cost
                BestSol = pop(i);
            end
        end

        %% Convergence
        convergenceCurve(t) = BestSol.Cost;

    end

    bestFitness = BestSol.Cost;
    bestPosition = BestSol.Position;

end

%% Roulette Wheel Selection (Local Function)
function i = RouletteWheelSelection(P)
    r = rand;
    C = cumsum(P);
    i = find(r <= C, 1, 'first');
end
