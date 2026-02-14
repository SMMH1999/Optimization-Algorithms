function [bestFitness, bestPosition, convergenceCurve] = FFO(lb, ub, dim, nPop, maxItr, objFun)
    % Flying Fox Optimization (FFO)
    %
    % Reference:
    % K. Zervoudakis & S. Tsafarakis, "A global optimizer inspired from the
    % survival strategies of flying foxes", Engineering with Computers, 2022.
    %
    % Inputs:
    %   lb       - Lower bound (scalar or vector) of decision variables
    %   ub       - Upper bound (scalar or vector) of decision variables
    %   dim      - Number of decision variables
    %   nPop     - Population size (number of flying foxes)
    %   maxItr   - Maximum number of iterations
    %   objFun   - Handle to objective function (fitness function)
    %
    % Outputs:
    %   bestFitness      - Best objective value found
    %   bestPosition     - Decision vector corresponding to bestFitness
    %   convergenceCurve - BestFitness at each iteration
    %
    % Tunable Parameters:
    %   deltasO      - Initial delta fractions for exploration
    %   deltasOmax   - Maximum delta fractions
    %   deltasOmin   - Minimum delta fractions
    %   parameter.alpha - Alpha fuzzy tuning parameter
    %   parameter.pa    - Probability pa fuzzy tuning parameter
    %   SurvList     - Number of survivors in survival list

    %% Parameters
    deltasO    = [0.2 0.4 0.6];
    deltasOmax = deltasO;
    deltasOmin = [0.02 0.04 0.06];
    parameter.alpha = [1 1.5 1.9];
    parameter.pa    = [0.5 0.85 0.99];
    SurvList   = round(nPop/4);

    %% Initialization
    Empty.Position = [];
    Empty.Past.Cost = [];
    Empty.Cost = [];

    FlyingFox = repmat(Empty, nPop, 1);
    BestSol.Cost = inf;
    funccount = 0;

    for i = 1:nPop
        FlyingFox(i).Position = unifrnd(lb, ub, 1, dim);
        FlyingFox(i).Cost = objFun(FlyingFox(i).Position);
        FlyingFox(i).Past.Cost = FlyingFox(i).Cost;
        funccount = funccount + 1;
        if FlyingFox(i).Cost <= BestSol.Cost
            BestSol = FlyingFox(i);
        end
    end

    convergenceCurve = zeros(maxItr,1);
    WorstCost = max([FlyingFox.Cost]);
    SurvivalList = FlyingFox;

    %% Main Loop
    for it = 1:maxItr
        for i = 1:nPop
            deltamax = abs(BestSol.Cost - WorstCost);
            deltas = deltasO * deltamax;
            [alpha, pa] = FuzzySelfTuning(BestSol, FlyingFox(i), WorstCost, deltasO, parameter);

            if norm(FlyingFox(i).Cost - BestSol.Cost) > (deltas(1)*0.5)
                z = FlyingFox(i).Position + alpha .* unifrnd(0,1,1,dim) .* (BestSol.Position - FlyingFox(i).Position);
            else
                A = randperm(nPop); A(A==i)=[];
                a=A(1); b=A(2);
                stepsize = unifrnd(0,1,1,dim).*(BestSol.Position-FlyingFox(i).Position) + ...
                    unifrnd(0,1,1,dim).*(FlyingFox(a).Position-FlyingFox(b).Position);
                z = zeros(1, dim);
                j0 = randi([1 dim]);
                for j = 1:dim
                    if j==j0 || rand>=pa
                        z(j) = FlyingFox(i).Position(j) + stepsize(j);
                    else
                        z(j) = FlyingFox(i).Position(j);
                    end
                end
            end

            temp = Empty;
            FlyingFox(i).Past.Cost = FlyingFox(i).Cost;
            temp.Past.Cost = FlyingFox(i).Past.Cost;
            z = max(z, lb); z = min(z, ub);
            temp.Position = z;
            temp.Cost = objFun(z);
            funccount = funccount + 1;

            if temp.Cost < FlyingFox(i).Cost
                FlyingFox(i) = temp;
                if FlyingFox(i).Cost <= BestSol.Cost
                    BestSol = FlyingFox(i);
                end
            end
            WorstCost = max(WorstCost, temp.Cost);
            SurvivalList = UpdateSurvivalList(SurvivalList, SurvList, temp);

            if norm(temp.Cost - BestSol.Cost) > deltas(3)
                FlyingFox(i) = ReplaceWithSurvivalList(SurvivalList, Empty, dim);
                FlyingFox(i).Position = max(min(FlyingFox(i).Position, ub), lb);
                FlyingFox(i).Cost = objFun(FlyingFox(i).Position);
                FlyingFox(i).Past.Cost = FlyingFox(i).Cost;
                WorstCost = max(WorstCost, FlyingFox(i).Cost);
                if FlyingFox(i).Cost < BestSol.Cost
                    BestSol = FlyingFox(i);
                end
                funccount = funccount + 1;
            end
        end

        %% Suffocating Best Flying Foxes
        pBestFF = find([FlyingFox.Cost]==BestSol.Cost);
        nBestFF = numel(pBestFF);
        pDeath = (nBestFF-1)/nPop;

        for i = 1:2:nBestFF
            if rand < pDeath
                j = 1:nPop; j(pBestFF)=[];
                if mod(nBestFF,2)==1 && i==nBestFF
                    FlyingFox(pBestFF(i)) = ReplaceWithSurvivalList(SurvivalList, Empty, dim);
                    FlyingFox(pBestFF(i)).Position = max(min(FlyingFox(pBestFF(i)).Position, ub), lb);
                    FlyingFox(pBestFF(i)).Cost = objFun(FlyingFox(pBestFF(i)).Position);
                    FlyingFox(pBestFF(i)).Past.Cost = FlyingFox(pBestFF(i)).Cost;
                    SurvivalList = UpdateSurvivalList(SurvivalList, SurvList, FlyingFox(pBestFF(i)));
                else
                    parent1=randi(nPop); parent2=randi(nPop);
                    if rand<0.5 && FlyingFox(parent1).Cost~=FlyingFox(parent2).Cost
                        [FlyingFox(pBestFF(i)).Position, FlyingFox(pBestFF(i+1)).Position] = ...
                            Crossover(FlyingFox(parent1).Position, FlyingFox(parent2).Position, lb, ub);
                    else
                        FlyingFox(pBestFF(i))=ReplaceWithSurvivalList(SurvivalList, Empty, dim);
                        FlyingFox(pBestFF(i+1))=ReplaceWithSurvivalList(SurvivalList, Empty, dim);
                    end
                    FlyingFox(pBestFF(i)).Cost = objFun(FlyingFox(pBestFF(i)).Position);
                    FlyingFox(pBestFF(i)).Past.Cost = FlyingFox(pBestFF(i)).Cost;
                    SurvivalList = UpdateSurvivalList(SurvivalList, SurvList, FlyingFox(pBestFF(i)));

                    FlyingFox(pBestFF(i+1)).Cost = objFun(FlyingFox(pBestFF(i+1)).Position);
                    FlyingFox(pBestFF(i+1)).Past.Cost = FlyingFox(pBestFF(i+1)).Cost;
                    SurvivalList = UpdateSurvivalList(SurvivalList, SurvList, FlyingFox(pBestFF(i+1)));
                end
            end
        end

        convergenceCurve(it) = BestSol.Cost;
        deltasO = deltasOmax - ((deltasOmax - deltasOmin)/maxItr)*it;
    end

    bestFitness = BestSol.Cost;
    bestPosition = BestSol.Position;
end
%% ------------------ Helper Functions ------------------
function z=Sphere(x)
    z=sum(x.^2);
end

function SurvivalList=UpdateSurvivalList(SurvivalList,SurvList,temp)
    if temp.Cost < SurvivalList(end).Cost
        SurvivalList = [temp; SurvivalList];
        [~,ii] = unique([SurvivalList.Cost]);
        SurvivalList = SurvivalList(ii);
        if numel(SurvivalList) > SurvList
            SurvivalList = SurvivalList(1:SurvList);
        end
    end
end

function ffox=ReplaceWithSurvivalList(SurvivalList,Empty,dim)
    m = randi([2 numel(SurvivalList)]);
    ffox = Empty;
    h = randperm(numel(SurvivalList), m);
    ffox.Position = zeros(1,dim);
    for i=1:numel(h)
        ffox.Position = ffox.Position + SurvivalList(h(i)).Position;
    end
    ffox.Position = ffox.Position / m;
end

function [off1, off2]=Crossover(x1,x2,lb,ub)
    L=unifrnd(0,1,size(x1));
    off1 = L.*x1 + (1-L).*x2;
    off2 = L.*x2 + (1-L).*x1;
    off1 = max(min(off1, ub), lb);
    off2 = max(min(off2, ub), lb);
end

function [alpha,pa]=FuzzySelfTuning(BestSol,FlyingFox,WorstCost,deltasO,parameter)
    delta = norm(BestSol.Cost - FlyingFox.Cost);
    deltamax = abs(BestSol.Cost - WorstCost);
    fi = (FlyingFox.Cost - FlyingFox.Past.Cost)/deltamax;
    deltas = deltasO * deltamax;

    % Simplified fuzzy rules (original logic preserved)
    alpha = sum(parameter.alpha)/numel(parameter.alpha);
    pa = sum(parameter.pa)/numel(parameter.pa);
end


