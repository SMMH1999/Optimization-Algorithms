function [bestFitness, bestPosition, convergenceCurve] = ALO(lb, ub, dim, nPop, maxItr, objFun)
    % Ant Lion Optimizer (ALO)
    %% Boundary handling
    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end
    lb = lb(:)';
    ub = ub(:)';

    %% Initialization
    antlionPos = rand(nPop, dim) .* (ub - lb) + lb;
    antPos     = rand(nPop, dim) .* (ub - lb) + lb;

    antlionFit = zeros(nPop,1);
    antFit     = zeros(nPop,1);

    for i = 1:nPop
        antlionFit(i) = objFun(antlionPos(i,:));
    end

    [antlionFit, idx] = sort(antlionFit);
    antlionPos = antlionPos(idx,:);

    elitePos = antlionPos(1,:);
    eliteFit = antlionFit(1);

    convergenceCurve = zeros(1,maxItr);
    convergenceCurve(1) = eliteFit;

    %% Main loop
    for t = 2:maxItr

        for i = 1:nPop

            % Roulette Wheel Selection (based on inverse fitness)
            weights = 1 ./ (antlionFit + eps);
            accumulation = cumsum(weights);
            p = rand * accumulation(end);
            index = find(accumulation >= p, 1, 'first');
            if isempty(index)
                index = 1;
            end

            % Random walks
            RA = randomWalk(antlionPos(index,:), t);
            RE = randomWalk(elitePos, t);

            antPos(i,:) = (RA + RE) / 2;
        end

        % Boundary control and fitness evaluation
        antPos = min(max(antPos, lb), ub);
        for i = 1:nPop
            antFit(i) = objFun(antPos(i,:));
        end

        % Combine populations
        doublePop = [antlionPos; antPos];
        doubleFit = [antlionFit; antFit];

        [doubleFit, idx] = sort(doubleFit);
        doublePop = doublePop(idx,:);

        antlionPos = doublePop(1:nPop,:);
        antlionFit = doubleFit(1:nPop);

        % Update elite
        if antlionFit(1) < eliteFit
            eliteFit = antlionFit(1);
            elitePos = antlionPos(1,:);
        end

        % Keep elite
        antlionPos(1,:) = elitePos;
        antlionFit(1)   = eliteFit;

        convergenceCurve(t) = eliteFit;
    end

    bestFitness  = eliteFit;
    bestPosition = elitePos;

    %% --- Local Random Walk Function ---
    function RW = randomWalk(antlion, currentIter)

        I = 1;
        if currentIter > maxItr/10
            I = 1 + 100*(currentIter/maxItr);
        end
        if currentIter > maxItr/2
            I = 1 + 1000*(currentIter/maxItr);
        end
        if currentIter > 3*maxItr/4
            I = 1 + 10000*(currentIter/maxItr);
        end
        if currentIter > 0.9*maxItr
            I = 1 + 100000*(currentIter/maxItr);
        end
        if currentIter > 0.95*maxItr
            I = 1 + 1000000*(currentIter/maxItr);
        end

        lb_t = lb ./ I;
        ub_t = ub ./ I;

        if rand < 0.5
            lb_t = lb_t + antlion;
        else
            lb_t = -lb_t + antlion;
        end

        if rand >= 0.5
            ub_t = ub_t + antlion;
        else
            ub_t = -ub_t + antlion;
        end

        RW = zeros(1,dim);

        for d = 1:dim
            X = [0 cumsum(2*(rand(maxItr,1)>0.5)-1)'];
            a = min(X);
            b = max(X);
            c = lb_t(d);
            d_ = ub_t(d);
            Xn = ((X - a) .* (d_ - c)) ./ (b - a + eps) + c;
            RW(d) = Xn(currentIter);
        end
    end

end
