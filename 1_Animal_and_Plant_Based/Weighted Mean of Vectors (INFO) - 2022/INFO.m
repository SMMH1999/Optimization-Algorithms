function [bestFitness, bestPosition, convergenceCurve] = INFO(lb, ub, dim, nPop, maxItr, objFun)
    %% INFO: Weighted Mean of Vectors Optimization Algorithm
    % INFO: An Efficient Optimization Algorithm based on Weighted Mean of Vectors
    % Authors: Iman Ahmadianfar, Ali Asghar Heidari, Saeed Noushadian, Huiling Chen, Amir H. Gandomi
    %
    % INPUTS:
    %   lb        : lower bound (scalar or vector) [1 x dim]
    %   ub        : upper bound (scalar or vector) [1 x dim]
    %   dim       : problem dimension
    %   nPop      : population size
    %   maxItr    : maximum number of iterations
    %   objFun    : handle of objective function, f(x)
    %
    % OUTPUTS:
    %   bestFitness       : best objective value found
    %   bestPosition      : best solution vector found
    %   convergenceCurve  : vector of bestFitness at each iteration

    %% Initialization
    Cost = zeros(nPop,1);
    M = zeros(nPop,1);

    X = initialization(nPop, dim, ub, lb);

    for i = 1:nPop
        Cost(i) = objFun(X(i,:));
        M(i) = Cost(i);
    end

    [~, ind] = sort(Cost);
    bestPosition = X(ind(1),:);
    bestFitness = Cost(ind(1));

    Worst_Cost = Cost(ind(end));
    Worst_X = X(ind(end),:);

    I = randi([2 5]);
    Better_X = X(ind(I),:);
    Better_Cost = Cost(ind(I));

    convergenceCurve = zeros(1,maxItr);

    %% Main Loop of INFO
    for it = 1:maxItr
        alpha = 2*exp(-4*(it/maxItr));

        M_Best = bestFitness;
        M_Better = Better_Cost;
        M_Worst = Worst_Cost;

        for i = 1:nPop
            % Updating rule stage
            del = 2*rand*alpha - alpha;
            sigm = 2*rand*alpha - alpha;

            % Select three random solutions
            A1 = randperm(nPop);
            A1(A1==i) = [];
            a = A1(1); b = A1(2); c = A1(3);

            e = 1e-25;
            epsi = e*rand;

            omg = max([M(a) M(b) M(c)]);
            MM = [(M(a)-M(b)) (M(a)-M(c)) (M(b)-M(c))];

            W(1) = cos(MM(1)+pi)*exp(-MM(1)/omg);
            W(2) = cos(MM(2)+pi)*exp(-MM(2)/omg);
            W(3) = cos(MM(3)+pi)*exp(-MM(3)/omg);
            Wt = sum(W);

            WM1 = del.*(W(1).*(X(a,:)-X(b,:)) + W(2).*(X(a,:)-X(c,:)) + W(3).*(X(b,:)-X(c,:)))/(Wt+1) + epsi;

            omg = max([M_Best M_Better M_Worst]);
            MM = [(M_Best-M_Better) (M_Best-M_Better) (M_Better-M_Worst)];

            W(1) = cos(MM(1)+pi)*exp(-MM(1)/omg);
            W(2) = cos(MM(2)+pi)*exp(-MM(2)/omg);
            W(3) = cos(MM(3)+pi)*exp(-MM(3)/omg);
            Wt = sum(W);

            WM2 = del.*(W(1).*(bestPosition-Better_X) + W(2).*(bestPosition-Worst_X) + W(3).*(Better_X-Worst_X))/(Wt+1) + epsi;

            % Determine MeanRule
            r = unifrnd(0.1,0.5);
            MeanRule = r.*WM1 + (1-r).*WM2;

            if rand < 0.5
                z1 = X(i,:) + sigm.*(rand.*MeanRule) + randn.*(bestPosition-X(a,:))./(M_Best-M(a)+1);
                z2 = bestPosition + sigm.*(rand.*MeanRule) + randn.*(X(a,:)-X(b,:))./(M(a)-M(b)+1);
            else
                z1 = X(a,:) + sigm.*(rand.*MeanRule) + randn.*(X(b,:)-X(c,:))./(M(b)-M(c)+1);
                z2 = Better_X + sigm.*(rand.*MeanRule) + randn.*(X(a,:)-X(b,:))./(M(a)-M(b)+1);
            end

            % Vector combining stage
            u = zeros(1,dim);
            for j = 1:dim
                mu = 0.05*randn;
                if rand < 0.5
                    if rand < 0.5
                        u(j) = z1(j) + mu*abs(z1(j)-z2(j));
                    else
                        u(j) = z2(j) + mu*abs(z1(j)-z2(j));
                    end
                else
                    u(j) = X(i,j);
                end
            end

            % Local search stage
            if rand < 0.5
                L = rand < 0.5;
                v1 = (1-L)*2*rand + L;
                v2 = rand.*L + (1-L);
                Xavg = (X(a,:) + X(b,:) + X(c,:))/3;
                phi = rand;
                Xrnd = phi.*Xavg + (1-phi).*(phi.*Better_X + (1-phi).*bestPosition);
                RandnVec = L.*randn(1,dim) + (1-L).*randn(1,dim);

                if rand < 0.5
                    u = bestPosition + RandnVec.*(MeanRule + randn.*(bestPosition-X(a,:)));
                else
                    u = Xrnd + RandnVec.*(MeanRule + randn.*(v1*bestPosition - v2*Xrnd));
                end
            end

            % Boundary check
            New_X = BC(u, lb, ub);
            New_Cost = objFun(New_X);

            if New_Cost < Cost(i)
                X(i,:) = New_X;
                Cost(i) = New_Cost;
                M(i) = New_Cost;
                if New_Cost < bestFitness
                    bestPosition = New_X;
                    bestFitness = New_Cost;
                end
            end
        end

        % Update worst and better solutions
        [~, ind] = sort(Cost);
        Worst_X = X(ind(end),:);
        Worst_Cost = Cost(ind(end));
        I = randi([2 5]);
        Better_X = X(ind(I),:);
        Better_Cost = Cost(ind(I));

        % Update convergence curve
        convergenceCurve(it) = bestFitness;
    end

end

%% Boundary control
function X = BC(X, lb, ub)
    Flag4ub = X > ub;
    Flag4lb = X < lb;
    X = X.*(~(Flag4ub + Flag4lb)) + ub.*Flag4ub + lb.*Flag4lb;
end

%% Initialization function
function X = initialization(nP, dim, ub, lb)
    if isscalar(ub)
        X = rand(nP, dim).*(ub-lb) + lb;
    else
        X = zeros(nP, dim);
        for i = 1:dim
            X(:,i) = rand(nP,1).*(ub(i)-lb(i)) + lb(i);
        end
    end
end
