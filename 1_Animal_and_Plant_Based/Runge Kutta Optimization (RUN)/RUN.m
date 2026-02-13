function [bestFitness, bestPosition, convergenceCurve] = RUN(lb, ub, dim, nPop, maxItr, objFun)
    %% RUN: Runge Kutta Optimizer
    % RUN Beyond the Metaphor: An Efficient Optimization Algorithm Based on Runge Kutta Method
    %
    % Author: Iman Ahmadianfar, Ali Asghar Heidari, Amir H. Gandomi, Xuefeng Chu, Huiling Chen
    %
    % INPUTS:
    %   lb       : Lower bounds (scalar or 1 x dim vector)
    %   ub       : Upper bounds (scalar or 1 x dim vector)
    %   dim      : Number of decision variables
    %   nPop     : Population size
    %   maxItr   : Maximum number of iterations
    %   objFun   : Handle to the objective function
    %
    % OUTPUTS:
    %   bestFitness       : Best objective function value found
    %   bestPosition      : Position vector of the best solution
    %   convergenceCurve  : Best fitness value at each iteration

    %% Initialization
    Cost = zeros(nPop,1);
    X = initializePopulation(nPop, dim, lb, ub);
    convergenceCurve = zeros(1, maxItr);

    for i = 1:nPop
        Cost(i) = objFun(X(i,:));
    end

    [bestFitness, ind] = min(Cost);
    bestPosition = X(ind,:);
    convergenceCurve(1) = bestFitness;

    %% Main Loop
    for it = 2:maxItr
        f = 20 * exp(-12 * (it/maxItr));
        Xavg = mean(X, 1);
        SF = 2 * (0.5 - rand(1, nPop)) * f;

        for i = 1:nPop
            [~, ind_l] = min(Cost);
            lBest = X(ind_l,:);

            [A, B, C] = getRandomIndices(nPop, i);
            [~, ind1] = min(Cost([A B C]));

            % Delta X
            gama = rand .* (X(i,:) - rand(1,dim) .* (ub - lb)) .* exp(-4*it/maxItr);
            Stp = rand(1,dim) .* ((bestPosition - rand .* Xavg) + gama);
            DelX = 2 * rand(1,dim) .* abs(Stp);

            % Determine Xb and Xw
            if Cost(i) < Cost(ind1)
                Xb = X(i,:);
                Xw = X(ind1,:);
            else
                Xb = X(ind1,:);
                Xw = X(i,:);
            end

            SM = rungeKutta(Xb, Xw, DelX);

            L = rand(1,dim) < 0.5;
            Xc = L .* X(i,:) + (1-L) .* X(A,:);
            Xm = L .* bestPosition + (1-L) .* lBest;

            vec = [1, -1];
            flag = floor(2*rand(1,dim)+1);
            r = vec(flag);

            g = 2*rand;
            mu = 0.5 + 0.1*randn(1,dim);

            if rand < 0.5
                Xnew = (Xc + r .* SF(i) .* g .* Xc) + SF(i) .* SM + mu .* (Xm - Xc);
            else
                Xnew = (Xm + r .* SF(i) .* g .* Xm) + SF(i) .* SM + mu .* (X(A,:) - X(B,:));
            end

            % Boundary control
            Xnew = min(max(Xnew, lb), ub);
            CostNew = objFun(Xnew);

            if CostNew < Cost(i)
                X(i,:) = Xnew;
                Cost(i) = CostNew;
            end

            %% Enhanced Solution Quality (ESQ)
            if rand < 0.5
                EXP = exp(-5*rand*it/maxItr);
                r_int = randi([-1,1]);

                u = 2*rand(1,dim);
                w = unifrnd(0,2,1,dim) .* EXP;

                [A,B,C] = getRandomIndices(nPop,i);
                XavgESQ = (X(A,:) + X(B,:) + X(C,:)) / 3;

                beta = rand(1,dim);
                Xnew1 = beta .* bestPosition + (1-beta) .* XavgESQ;

                Xnew2 = zeros(1,dim);
                for j=1:dim
                    if w(j) < 1
                        Xnew2(j) = Xnew1(j) + r_int*w(j)*abs((Xnew1(j)-XavgESQ(j)) + randn);
                    else
                        Xnew2(j) = (Xnew1(j) - XavgESQ(j)) + r_int*w(j)*abs((u(j)*Xnew1(j)-XavgESQ(j)) + randn);
                    end
                end

                Xnew2 = min(max(Xnew2, lb), ub);
                CostNew = objFun(Xnew2);

                if CostNew < Cost(i)
                    X(i,:) = Xnew2;
                    Cost(i) = CostNew;
                else
                    if rand < w(randi(dim))
                        SM = rungeKutta(X(i,:), Xnew2, DelX);
                        Xnew = (Xnew2 - rand .* Xnew2) + SF(i) * (SM + (2*rand(1,dim).*bestPosition - Xnew2));
                        Xnew = min(max(Xnew, lb), ub);
                        CostNew = objFun(Xnew);
                        if CostNew < Cost(i)
                            X(i,:) = Xnew;
                            Cost(i) = CostNew;
                        end
                    end
                end
            end

            % Update global best
            if Cost(i) < bestFitness
                bestPosition = X(i,:);
                bestFitness = Cost(i);
            end
        end

        convergenceCurve(it) = bestFitness;
    end

end

%% Supporting Functions

function X = initializePopulation(nP, dim, lb, ub)
    X = rand(nP, dim) .* (ub - lb) + lb;
end

function [A,B,C] = getRandomIndices(nP, i)
    Qi = randperm(nP); Qi(Qi==i) = [];
    A = Qi(1); B = Qi(2); C = Qi(3);
end

function SM = rungeKutta(XB, XW, DelX)
    dim = numel(XB);
    C = randi([1 2])*(1-rand);
    r1 = rand(1,dim); r2 = rand(1,dim);

    K1 = 0.5*(rand*XW - C.*XB);
    K2 = 0.5*(rand*(XW + r2.*K1.*DelX/2) - (C*XB + r1.*K1.*DelX/2));
    K3 = 0.5*(rand*(XW + r2.*K2.*DelX/2) - (C*XB + r1.*K2.*DelX/2));
    K4 = 0.5*(rand*(XW + r2.*K3.*DelX) - (C*XB + r1.*K3.*DelX));

    SM = (K1 + 2*K2 + 2*K3 + K4)/6;
end

function z = unifrnd(a,b,c,dim)
    a2 = a/2; b2 = b/2;
    mu = a2 + b2; sig = b2 - a2;
    z = mu + sig .* (2*rand(c,dim) - 1);
end
