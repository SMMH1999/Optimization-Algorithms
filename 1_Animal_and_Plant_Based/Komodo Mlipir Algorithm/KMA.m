function [bestFitness, bestPosition, convergenceCurve] = KMA(lb, ub, dim, nPop, maxItr, objFun)
    %==========================================================================
    % Komodo Mlipir Algorithm (KMA) version 1.0 - Unified MATLAB function
    %
    % Developed by S. Suyanto et al., Applied Soft Computing, 2021
    %
    % DESCRIPTION:
    %   Population-based metaheuristic inspired by Komodo dragon hunting behavior.
    %   Performs global optimization with adaptive population size and gender roles.
    %
    % INPUTS:
    %   lb        - lower bounds (1xD vector or scalar)
    %   ub        - upper bounds (1xD vector or scalar)
    %   dim       - number of decision variables
    %   nPop      - initial population size
    %   maxItr    - maximum number of iterations
    %   objFun    - handle to objective function: f = objFun(x)
    %
    % OUTPUTS:
    %   bestFitness      - best objective function value found
    %   bestPosition     - corresponding decision variable vector
    %   convergenceCurve - best fitness at each iteration
    %==========================================================================

    %% -------------------- Initialization -----------------------------------
    if isscalar(lb), lb = repmat(lb,1,dim); end
    if isscalar(ub), ub = repmat(ub,1,dim); end

    Pop = repmat(lb,nPop,1) + rand(nPop,dim) .* (repmat(ub-lb,nPop,1));  % Initial population
    FX = zeros(nPop,1);
    for i = 1:nPop
        FX(i) = objFun(Pop(i,:));
    end
    [FX, idx] = sort(FX);
    Pop = Pop(idx,:);
    bestPosition = Pop(1,:);
    bestFitness  = FX(1);

    NumBM      = floor(nPop/2);     % Big males
    MlipirRate = (dim-1)/dim;       % Small males movement rate
    MutRate    = 0.5;               % Female mutation probability
    MutRadius  = 0.5;               % Female mutation step factor
    convergenceCurve = zeros(1,maxItr);

    %% -------------------- Main Loop -----------------------------------------
    for iter = 1:maxItr
        % Partition population
        BigMales    = Pop(1:NumBM,:);
        BigMalesFX  = FX(1:NumBM);
        Female      = Pop(NumBM+1,:);
        FemaleFX    = FX(NumBM+1);
        SmallMales  = Pop(NumBM+2:end,:);
        SmallMalesFX= FX(NumBM+2:end);

        % --- Move Big Males and Female ---
        [BigMales, BigMalesFX, Female, FemaleFX] = MoveBigMalesFemale(BigMales, BigMalesFX, Female, FemaleFX, dim, lb, ub, objFun, MutRate, MutRadius);

        % --- Move Small Males (Mlipir) ---
        [SmallMales, SmallMalesFX] = MoveSmallMales(SmallMales, SmallMalesFX, BigMales, MlipirRate, dim, lb, ub, objFun);

        % --- Update population and sort ---
        Pop = [BigMales; Female; SmallMales];
        FX  = [BigMalesFX; FemaleFX; SmallMalesFX];
        [FX, idx] = sort(FX);
        Pop = Pop(idx,:);

        % --- Track best ---
        bestPosition = Pop(1,:);
        bestFitness  = FX(1);
        convergenceCurve(iter) = bestFitness;
    end
end

%% -------------------- Helper Functions ----------------------------------

function [BigMales, BigMalesFX, Female, FemaleFX] = MoveBigMalesFemale(BigMales, BigMalesFX, Female, FemaleFX, dim, lb, ub, objFun, MutRate, MutRadius)
    TempBM   = BigMales;
    TempBMFX = BigMalesFX;
    for i = 1:size(BigMales,1)
        VM = zeros(1,dim);
        idxs = randperm(size(BigMales,1));
        maxFollow = randi(2);
        fol = 0;
        for j = 1:length(idxs)
            if idxs(j) ~= i
                if BigMalesFX(idxs(j)) < TempBMFX(i) || rand<0.5
                    VM = VM + rand*(BigMales(idxs(j),:) - TempBM(i,:));
                else
                    VM = VM + rand*(TempBM(i,:) - BigMales(idxs(j),:));
                end
            end
            fol = fol + 1;
            if fol >= maxFollow, break; end
        end
        NewBM = TempBM(i,:) + VM;
        NewBM = max(min(NewBM, ub), lb);
        TempBM(i,:) = NewBM;
        TempBMFX(i) = objFun(NewBM);
    end
    [BigMales, BigMalesFX] = Replacement(BigMales, BigMalesFX, TempBM, TempBMFX);

    % --- Female reproduction ---
    WinnerBM = BigMales(1,:);
    WinnerFX = BigMalesFX(1);
    if WinnerFX < FemaleFX || rand<0.5
        Offsprings = CrossOver(WinnerBM, Female);
        fx1 = objFun(Offsprings(1,:));
        fx2 = objFun(Offsprings(2,:));
        if fx1 < fx2
            if fx1 < FemaleFX
                Female = Offsprings(1,:);
                FemaleFX = fx1;
            end
        else
            if fx2 < FemaleFX
                Female = Offsprings(2,:);
                FemaleFX = fx2;
            end
        end
    else
        NewFemale = Mutation(Female, dim, lb, ub, MutRate, MutRadius, objFun);
        if NewFemale.fx < FemaleFX
            Female = NewFemale.x;
            FemaleFX = NewFemale.fx;
        end
    end
end

function [SmallMales, SmallMalesFX] = MoveSmallMales(SmallMales, SmallMalesFX, BigMales, MlipirRate, dim, lb, ub, objFun)
    TempSM = SmallMales;
    TempSMFX = SmallMalesFX;
    HQ = BigMales;
    for i = 1:size(SmallMales,1)
        VM = zeros(1,dim);
        idxs = randperm(size(HQ,1));
        maxFollow = randi(1);
        fol = 0;
        for j = 1:length(idxs)
            D = max(1, min(dim-1, round(MlipirRate*dim)));
            moveIdx = idxs(1:D);
            B = zeros(1,dim); B(moveIdx) = 1;
            VM = VM + rand(1,dim).*(HQ(idxs(j),:).*B - SmallMales(i,:).*B);
            fol = fol + 1;
            if fol >= maxFollow, break; end
        end
        NewSM = SmallMales(i,:) + VM;
        TempSM(i,:) = max(min(NewSM, ub), lb);
        TempSMFX(i) = objFun(TempSM(i,:));
    end
    SmallMales = TempSM;
    SmallMalesFX = TempSMFX;
end

function [Z, FZ] = Replacement(X, FX, Y, FY)
    XY = [X;Y];
    FXFY = [FX;FY];
    [FXsorted, idx] = sort(FXFY);
    Z = XY(idx(1:size(X,1)),:);
    FZ = FXsorted(1:size(X,1));
end

function Offsprings = CrossOver(P1,P2)
    dim = length(P1);
    Offsprings = zeros(2,dim);
    for i=1:dim
        rval = rand;
        Offsprings(1,i) = rval*P1(i) + (1-rval)*P2(i);
        Offsprings(2,i) = rval*P2(i) + (1-rval)*P1(i);
    end
end

function NewFemale = Mutation(Female, dim, lb, ub, MutRate, MutRadius, objFun)
    MaxStep = MutRadius*(ub-lb);
    NewX = Female;
    for i=1:dim
        if rand<MutRate
            NewX(i) = NewX(i) + (2*rand-1)*MaxStep(i);
        end
    end
    NewX = max(min(NewX, ub), lb);
    NewFemale.x = NewX;
    NewFemale.fx = objFun(NewX);
end
