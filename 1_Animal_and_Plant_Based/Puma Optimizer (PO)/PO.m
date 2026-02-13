function [bestFitness, bestPosition, convergenceCurve] = PO(lb, ub, dim, nPop, maxItr, objFun)
    %_______________________________________________________________________________________________
    % Puma Optimizer Algorithm (POA)
    %
    % Source: Cluster Computing, DOI: 10.1007/s10586-023-04221-5
    % Authors: B. Abdollahzadeh, N. Khodadadi, S. Barshandeh, P. Trojovský,
    %          F. Soleimanian Gharehchopogh, E.M. El-kenawy, L. Abualigah, S. Mirjalili
    %
    % Developed in MATLAB R2011b-R2021a
    %
    % INPUTS:
    %   lb        : Lower bound (scalar or vector, 1 x dim)
    %   ub        : Upper bound (scalar or vector, 1 x dim)
    %   dim       : Dimension of the problem
    %   nPop      : Number of candidate solutions (population)
    %   maxItr    : Maximum number of iterations
    %   objFun    : Handle to the objective function: Cost = objFun(X)
    %
    % OUTPUTS:
    %   bestFitness     : Best cost obtained
    %   bestPosition    : Best position vector (1 x dim)
    %   convergenceCurve: Vector of best cost per iteration
    %
    %_______________________________________________________________________________________________

    %% Parameter Initialization
    UnSelected = ones(1,2); % 1: Exploration, 2: Exploitation
    F3_Explore = 0;
    F3_Exploit = 0;
    Seq_Time_Explore = ones(1,3);
    Seq_Time_Exploit = ones(1,3);
    Seq_Cost_Explore = ones(1,3);
    Seq_Cost_Exploit = ones(1,3);
    Score_Explore = 0;
    Score_Exploit = 0;
    PF = [0.5 0.5 0.3];  % Intensification (F1,F2) & Diversification (F3)
    PF_F3 = [];
    Mega_Explor = 0.99;
    Mega_Exploit = 0.99;

    %% Initialize Population
    Sol(nPop) = struct('X', [], 'Cost', []);
    for i = 1:nPop
        Sol(i).X = unifrnd(lb, ub, 1, dim);
        Sol(i).Cost = objFun(Sol(i).X);
    end
    [~, ind] = min([Sol.Cost]);
    Best = Sol(ind);
    Initial_Best = Best;
    Flag_Change = 1;
    convergenceCurve = zeros(1, maxItr);

    %% Unexperienced Phase (first 3 iterations)
    Costs_Explor = zeros(1,3);
    Costs_Exploit = zeros(1,3);
    for Iter = 1:3
        Sol_Explor = Exploration(Sol, lb, ub, dim, nPop, objFun);
        Costs_Explor(Iter) = min([Sol_Explor.Cost]);

        Sol_Exploit = Exploitation(Sol, lb, ub, dim, nPop, Best, maxItr, Iter, objFun);
        Costs_Exploit(Iter) = min([Sol_Exploit.Cost]);

        Sol = [Sol Sol_Explor Sol_Exploit];
        [~, sind] = sort([Sol.Cost]);
        Sol = Sol(sind(1:nPop));
        Best = Sol(1);

        convergenceCurve(Iter) = Best.Cost;
    end

    %% Hyper Initialization
    Seq_Cost_Explore(1) = abs(Initial_Best.Cost - Costs_Explor(1));
    Seq_Cost_Exploit(1) = abs(Initial_Best.Cost - Costs_Exploit(1));
    Seq_Cost_Explore(2) = abs(Costs_Explor(2) - Costs_Explor(1));
    Seq_Cost_Exploit(2) = abs(Costs_Exploit(2) - Costs_Exploit(1));
    Seq_Cost_Explore(3) = abs(Costs_Explor(3) - Costs_Explor(2));
    Seq_Cost_Exploit(3) = abs(Costs_Exploit(3) - Costs_Exploit(2));

    for i = 1:3
        if Seq_Cost_Explore(i) ~= 0
            PF_F3 = [PF_F3, Seq_Cost_Explore(i)];
        end
        if Seq_Cost_Exploit(i) ~= 0
            PF_F3 = [PF_F3, Seq_Cost_Exploit(i)];
        end
    end

    % Compute initial scores
    F1_Explor = PF(1)*(Seq_Cost_Explore(1)/Seq_Time_Explore(1));
    F1_Exploit = PF(1)*(Seq_Cost_Exploit(1)/Seq_Time_Exploit(1));
    F2_Explor = PF(2)*sum(Seq_Cost_Explore)/sum(Seq_Time_Explore);
    F2_Exploit = PF(2)*sum(Seq_Cost_Exploit)/sum(Seq_Time_Exploit);
    Score_Explore = F1_Explor + F2_Explor;
    Score_Exploit = F1_Exploit + F2_Exploit;

    %% Experienced Phase (Iterations 4:maxItr)
    for Iter = 4:maxItr
        if Score_Explore > Score_Exploit
            SelectFlag = 1;
            Sol = Exploration(Sol, lb, ub, dim, nPop, objFun);
            Count_select = UnSelected;
            UnSelected(2) = UnSelected(2) + 1;
            UnSelected(1) = 1;
            F3_Explore = PF(3);
            F3_Exploit = F3_Exploit + PF(3);
            [~, TBind] = min([Sol.Cost]);
            TBest = Sol(TBind);
            Seq_Cost_Explore = [abs(Best.Cost - TBest.Cost), Seq_Cost_Explore(1:2)];
            if Seq_Cost_Explore(1) ~= 0
                PF_F3 = [PF_F3, Seq_Cost_Explore(1)];
            end
            if TBest.Cost < Best.Cost
                Best = TBest;
            end
        else
            SelectFlag = 2;
            Sol = Exploitation(Sol, lb, ub, dim, nPop, Best, maxItr, Iter, objFun);
            Count_select = UnSelected;
            UnSelected(1) = UnSelected(1) + 1;
            UnSelected(2) = 1;
            F3_Explore = F3_Explore + PF(3);
            F3_Exploit = PF(3);
            [~, TBind] = min([Sol.Cost]);
            TBest = Sol(TBind);
            Seq_Cost_Exploit = [abs(Best.Cost - TBest.Cost), Seq_Cost_Exploit(1:2)];
            if Seq_Cost_Exploit(1) ~= 0
                PF_F3 = [PF_F3, Seq_Cost_Exploit(1)];
            end
            if TBest.Cost < Best.Cost
                Best = TBest;
            end
        end

        if Flag_Change ~= SelectFlag
            Flag_Change = SelectFlag;
            Seq_Time_Explore = [Count_select(1), Seq_Time_Explore(1:2)];
            Seq_Time_Exploit = [Count_select(2), Seq_Time_Exploit(1:2)];
        end

        % Hyper Initialization
        F1_Explor = PF(1)*(Seq_Cost_Explore(1)/Seq_Time_Explore(1));
        F1_Exploit = PF(1)*(Seq_Cost_Exploit(1)/Seq_Time_Exploit(1));
        F2_Explor = PF(2)*sum(Seq_Cost_Explore)/sum(Seq_Time_Explore);
        F2_Exploit = PF(2)*sum(Seq_Cost_Exploit)/sum(Seq_Time_Exploit);

        % Mega adjustment
        if Score_Explore < Score_Exploit
            Mega_Explor = max(Mega_Explor-0.01, 0.01);
            Mega_Exploit = 0.99;
        elseif Score_Explore > Score_Exploit
            Mega_Explor = 0.99;
            Mega_Exploit = max(Mega_Exploit-0.01, 0.01);
        end
        lmn_Explore = 1 - Mega_Explor;
        lmn_Exploit = 1 - Mega_Exploit;

        Score_Explore = Mega_Explor*(F1_Explor+F2_Explor) + lmn_Explore*min(PF_F3)*F3_Explore;
        Score_Exploit = Mega_Exploit*(F1_Exploit+F2_Exploit) + lmn_Exploit*min(PF_F3)*F3_Exploit;

        convergenceCurve(Iter) = Best.Cost;
    end

    bestFitness = Best.Cost;
    bestPosition = Best.X;

    %% ---------------- Nested Functions ----------------
    function SolOut = Exploration(SolIn, lb, ub, dim, nSol, objFun)
        pCR = 0.2;
        PCR = 1 - pCR;
        p = PCR / nSol;
        [~, sind] = sort([SolIn.Cost]);
        SolIn = SolIn(sind);
        SolOut = SolIn;
        for i = 1:nSol
            x = SolIn(i).X;
            A = randperm(nSol);
            A(A == i) = [];
            a = A(1); b = A(2); c = A(3); d = A(4); e = A(5); f = A(6);
            G = 2*rand - 1;
            if rand < 0.5
                y = rand(1,dim).*(ub - lb) + lb;
            else
                y = SolIn(a).X + G.*(SolIn(a).X - SolIn(b).X) + G.*(((SolIn(a).X-SolIn(b).X)-(SolIn(c).X-SolIn(d).X))+((SolIn(c).X-SolIn(d).X)-(SolIn(e).X-SolIn(f).X)));
            end
            y = max(y, lb); y = min(y, ub);
            z = zeros(size(x));
            j0 = randi([1 numel(x)]);
            for j = 1:numel(x)
                if j == j0 || rand <= pCR
                    z(j) = y(j);
                else
                    z(j) = x(j);
                end
            end
            NewSol.X = z;
            NewSol.Cost = objFun(z);
            if NewSol.Cost < SolOut(i).Cost
                SolOut(i) = NewSol;
            else
                pCR = pCR + p;
            end
        end
    end

    function SolOut = Exploitation(SolIn, lb, ub, dim, nSol, BestSol, maxItr, Iter, objFun)
        Q = 0.67; Beta = 2;
        SolOut = SolIn;
        mbest = mean(reshape([SolIn.X], dim, nSol),2)';
        for i = 1:nSol
            beta1 = 2*rand;
            beta2 = randn(1,dim);
            w = randn(1,dim);
            v = randn(1,dim);
            F1 = randn(1,dim).*exp(2 - Iter*(2/maxItr));
            F2 = w.*v.^2.*cos(2*rand*w);
            R_1 = 2*rand-1;
            S1 = 2*rand-1 + randn(1,dim);
            S2 = F1*R_1.*SolIn(i).X + F2*(1-R_1).*BestSol.X;
            VEC = S2./S1;
            if rand() <= 0.5
                if rand() > Q
                    Xnew = BestSol.X + beta1.*(exp(beta2).*(SolIn(randi(nSol,1)).X - SolIn(i).X));
                else
                    Xnew = beta1.*VEC - BestSol.X;
                end
            else
                r1 = randi(nSol);
                Xnew = (mbest.*SolIn(r1).X - ((-1)^(randi([0 1],1,1))).*SolIn(i).X)./(1 + Beta*rand);
            end
            Xnew = max(min(Xnew, ub), lb);
            NewSol.X = Xnew;
            NewSol.Cost = objFun(Xnew);
            if NewSol.Cost < SolOut(i).Cost
                SolOut(i) = NewSol;
            end
        end
    end

end
