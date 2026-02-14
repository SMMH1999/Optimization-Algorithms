function [bestFitness, bestPosition, convergenceCurve] = KOA(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Kepler Optimization Algorithm (KOA)
    % =========================================================================
    % Author: Reda Mohamed & Mohamed Abdel-Basset
    % Source Paper:
    % Abdel-Basset, M., Mohamed, R.
    % "Kepler optimization algorithm: A new metaheuristic algorithm inspired by
    % Kepler’s laws of planetary motion"
    % Knowledge-Based Systems, 2023.
    % DOI: 10.1016/j.knosys.2023.110454
    %
    % =========================================================================
    % Function Signature:
    % [bestFitness, bestPosition, convergenceCurve] = KOA(lb, ub, dim, nPop, maxItr, objFun)
    %
    % =========================================================================
    % Inputs:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Number of decision variables
    % nPop      : Number of search agents (planets)
    % maxItr    : Maximum number of function evaluations
    % objFun    : Objective function handle (minimization)
    %
    % Outputs:
    % bestFitness      : Best objective value found
    % bestPosition     : Best solution vector found (1×dim)
    % convergenceCurve : Best fitness value at each iteration
    %
    % =========================================================================
    % Tunable Parameters:
    % Tc     : Cycle control parameter (default = 3)
    % M0     : Initial mass constant (default = 0.1)
    % lambda : Mass decay coefficient (default = 15)
    % =========================================================================

    %% --------------------- Parameter Handling -----------------------------
    if numel(lb) == 1
        lb = repmat(lb, 1, dim);
    end
    if numel(ub) == 1
        ub = repmat(ub, 1, dim);
    end

    %% --------------------- Initialization ---------------------------------
    Tc = 3;
    M0 = 0.1;
    lambda = 15;

    bestPosition = zeros(1, dim);
    bestFitness = inf;
    convergenceCurve = zeros(1, maxItr);

    orbital = rand(1, nPop);
    T = abs(randn(1, nPop));
    Positions = initialization(nPop, dim, ub, lb);

    fitness = zeros(1, nPop);

    for i = 1:nPop
        fitness(i) = objFun(Positions(i,:));
        if fitness(i) < bestFitness
            bestFitness = fitness(i);
            bestPosition = Positions(i,:);
        end
    end

    %% --------------------- Main Loop --------------------------------------
    t = 1;
    while t <= maxItr

        [sortedFit, idx] = sort(fitness);
        worstFitness = sortedFit(end);

        M = M0 * exp(-lambda * (t / maxItr));

        % Distance R
        R = sqrt(sum((Positions - bestPosition).^2, 2))';

        sumFit = sum(fitness - worstFitness);
        if sumFit == 0
            sumFit = eps;
        end

        MS = rand(1,nPop) .* (bestFitness - worstFitness) ./ sumFit;
        m = (fitness - worstFitness) ./ sumFit;

        % Normalize
        Rnorm = (R - min(R)) ./ (max(R) - min(R) + eps);
        MSnorm = (MS - min(MS)) ./ (max(MS) - min(MS) + eps);
        Mnorm = (m - min(m)) ./ (max(m) - min(m) + eps);

        Fg = orbital .* M .* ((MSnorm .* Mnorm) ./ (Rnorm.^2 + eps)) + rand(1,nPop);

        a1 = rand(1,nPop) .* (T.^2 .* (M .* (MS + m)) ./ (4*pi*pi)).^(1/3);

        for i = 1:nPop

            a2 = -1 + -1*(rem(t, maxItr/Tc)/(maxItr/Tc));
            n = (a2 - 1)*rand + 1;

            a = randi(nPop);
            b = randi(nPop);

            rd = rand(1,dim);
            r = rand;
            U1 = rd < r;
            O_P = Positions(i,:);

            if rand < rand
                h = 1/(exp(n*randn));
                Xm = (Positions(b,:) + bestPosition + Positions(i,:)) / 3;
                Positions(i,:) = Positions(i,:).*U1 + ...
                    (Xm + h.*(Xm - Positions(a,:))).*(~U1);
            else

                if rand < 0.5
                    f = 1;
                else
                    f = -1;
                end

                L = sqrt(M*(MS(i)+m(i))*abs((2/(R(i)+eps)) - (1/(a1(i)+eps))));
                U = rd > rand(1,dim);

                if Rnorm(i) < 0.5
                    Mtemp = rand*(1-r)+r;
                    l = L*Mtemp.*U;
                    Mv = rand*(1-rd)+rd;
                    l1 = L.*Mv.*(~U);

                    V = l.*(2*rand*Positions(i,:) - Positions(a,:)) + ...
                        l1.*(Positions(b,:) - Positions(a,:)) + ...
                        (1-Rnorm(i))*f.*U1.*rand(1,dim).*(ub-lb);
                else
                    U2 = rand > rand;
                    V = rand*L.*(Positions(a,:) - Positions(i,:)) + ...
                        (1-Rnorm(i))*f*U2.*rand(1,dim).*(rand*ub - lb);
                end

                if rand < 0.5
                    f = 1;
                else
                    f = -1;
                end

                Positions(i,:) = (Positions(i,:) + V.*f) + ...
                    (Fg(i) + abs(randn)).*U.*(bestPosition - Positions(i,:));
            end

            % Boundary control
            if rand < rand
                for j = 1:dim
                    if Positions(i,j) > ub(j) || Positions(i,j) < lb(j)
                        Positions(i,j) = lb(j) + rand*(ub(j)-lb(j));
                    end
                end
            else
                Positions(i,:) = min(max(Positions(i,:), lb), ub);
            end

            % Evaluation
            newFit = objFun(Positions(i,:));

            if newFit < fitness(i)
                fitness(i) = newFit;
                if newFit < bestFitness
                    bestFitness = newFit;
                    bestPosition = Positions(i,:);
                end
            else
                Positions(i,:) = O_P;
            end

            convergenceCurve(t) = bestFitness;

            t = t + 1;
            if t > maxItr
                break;
            end
        end
    end

end

%% ======================= Initialization Function ========================
function Positions = initialization(nPop, dim, ub, lb)

    Boundary_no = length(ub);

    if Boundary_no == 1
        Positions = rand(nPop, dim).*(ub - lb) + lb;
    else
        Positions = zeros(nPop, dim);
        for i = 1:dim
            Positions(:,i) = rand(nPop,1).*(ub(i) - lb(i)) + lb(i);
        end
    end

end
