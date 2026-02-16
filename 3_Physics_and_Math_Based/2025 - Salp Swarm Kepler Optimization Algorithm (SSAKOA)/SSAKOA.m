function [bestFitness, bestPosition, convergenceCurve] = SSAKOA(lb, ub, dim, nPop, maxItr, objFun)
    %__________________________________________________________________________________%
    % SSAKOA: Hybrid Salp Swarm Kepler Optimization Algorithm
    % Author: Dr. Aykut Fatih GÜVEN (afatih.guven@yalova.edu.tr)
    % Reference: A novel hybrid Salp Swarm Kepler optimization for optimal sizing and
    %            energy management of renewable microgrids with EV integration,
    %            Energy, 2025, 334:137696. DOI: 10.1016/j.energy.2025.137696
    %__________________________________________________________________________________%
    %
    % INPUTS:
    %   lb        : Lower bounds [1xdim or scalar]
    %   ub        : Upper bounds [1xdim or scalar]
    %   dim       : Dimension of the problem [scalar]
    %   nPop      : Population size [scalar]
    %   maxItr    : Maximum number of iterations [scalar]
    %   objFun    : Handle to the objective function [function handle]
    %
    % OUTPUTS:
    %   bestFitness       : Best fitness value found [scalar]
    %   bestPosition      : Best position vector found [1xdim]
    %   convergenceCurve  : Best fitness at each iteration [1xmaxItr]
    %
    % TUNABLE PARAMETERS:
    %   Tc     : Cycle controlling factor for KOA component
    %   M0     : Initial mass constant for KOA
    %   lambda : Decay factor for mass over iterations
    %__________________________________________________________________________________%

    %%---------------------- Handle bounds ----------------------%%
    if numel(lb) == 1
        lb = lb * ones(1,dim);
        ub = ub * ones(1,dim);
    elseif numel(lb) < dim
        lb = lb .* ones(1,dim);
        ub = ub .* ones(1,dim);
    end

    %%---------------------- Initialization --------------------%%
    convergenceCurve = zeros(1,maxItr);
    SalpPositions = initialization(nPop,dim,ub,lb);

    bestPosition = zeros(1,dim);
    bestFitness = inf;

    % KOA controlling parameters
    Tc = 3;
    M0 = 0.1;
    lambda = 15;

    SearchAgents_no = nPop;
    orbital = rand(1,SearchAgents_no);         % Orbital eccentricity
    T = abs(randn(1,SearchAgents_no));        % Orbital period
    Positions = SalpPositions;

    % Evaluate initial fitness
    SalpFitness = zeros(1,nPop);
    for i = 1:nPop
        SalpFitness(i) = objFun(SalpPositions(i,:));
    end

    [sortedSalpsFitness,sortedIndexes] = sort(SalpFitness);
    SortedSalps = SalpPositions(sortedIndexes,:);
    bestPosition = SortedSalps(1,:);
    bestFitness = sortedSalpsFitness(1);

    convergenceCurve(1) = bestFitness;

    %%------------------------ Main Loop -----------------------%%
    for l = 2:maxItr
        c1 = 2*exp(-(4*l/maxItr)^2); % Salp controlling factor

        %% KOA Component
        PL_Fit = SalpFitness;
        worstFitness = max(PL_Fit);
        M = M0 * exp(-lambda * (l-1)/maxItr);

        R = sqrt(sum((Positions - bestPosition).^2,2))'; % Distance to best
        sum_diff = sum(PL_Fit - worstFitness);
        MS = rand(1,SearchAgents_no) .* (bestFitness - worstFitness) / sum_diff;
        m  = (PL_Fit - worstFitness) / sum_diff;

        Rnorm  = (R - min(R)) / (max(R) - min(R) + eps);
        MSnorm = (MS - min(MS)) / (max(MS) - min(MS) + eps);
        Mnorm  = (m - min(m)) / (max(m) - min(m) + eps);

        Fg = orbital .* M .* ((MSnorm .* Mnorm) ./ (Rnorm.^2 + eps)) + rand(1,SearchAgents_no);
        a1 = rand(1,SearchAgents_no) .* (T.^2 .* (M .* (MS + m) / (4*pi^2))).^(1/3);

        %% Update Salp Positions
        for i = 1:SearchAgents_no
            a2 = -1 - 1 * (rem(l-1,maxItr/Tc)/(maxItr/Tc));
            n = (a2 - 1) * rand + 1;
            a_idx = randi(SearchAgents_no);
            b_idx = randi(SearchAgents_no);
            rd = rand(1,dim);
            r_scalar = rand;
            U1 = rd < r_scalar;

            if rand < rand
                h = 1 ./ exp(n .* randn(1,dim));
                Xm = (Positions(b_idx,:) + bestPosition + Positions(i,:))/3;
                Positions(i,:) = Positions(i,:) .* U1 + (Xm + h .* (Xm - Positions(a_idx,:))) .* (1 - U1);
            else
                f = 2*(rand<0.5)-1;
                L = sqrt(M*(MS(i)+m(i))*abs(2./(R(i)+eps) - 1./(a1(i)+eps)));
                U = rd > rand(1,dim);
                if Rnorm(i) < 0.5
                    lnew = L*M*U;
                    Mv = rand*(1-rd)+rd;
                    l1 = L.*Mv.*(1-U);
                    V = lnew.*(2*rand*Positions(i,:) - Positions(a_idx,:)) + l1.*(Positions(b_idx,:) - Positions(a_idx,:)) + (1-Rnorm(i))*f*U1.*rand(1,dim).*(ub-lb);
                else
                    U2 = rand > rand;
                    V = rand.*L.*(Positions(a_idx,:) - Positions(i,:)) + (1-Rnorm(i))*f*U2.*rand(1,dim).*(rand*ub - lb);
                end
                f = 2*(rand<0.5)-1;
                Positions(i,:) = (Positions(i,:) + V*f) + (Fg(i)+abs(randn)) * U .* (bestPosition - Positions(i,:));
            end

            %% Enforce boundaries
            Positions(i,:) = min(max(Positions(i,:), lb), ub);

            %% Salp update (SSA component)
            if i <= SearchAgents_no/2
                c2 = rand(1,dim); c3 = rand(1,dim);
                SalpPositions(i,:) = bestPosition + ((-1).^(c3<0.5)) .* c1 .* ((ub-lb).*c2 + lb);
            else
                SalpPositions(i,:) = Positions(i,:);
            end
        end

        %% Evaluate fitness and update best
        for i = 1:SearchAgents_no
            SalpPositions(i,:) = min(max(SalpPositions(i,:), lb), ub);
            fit_val = objFun(SalpPositions(i,:));
            if fit_val < SalpFitness(i)
                SalpFitness(i) = fit_val;
                if fit_val < bestFitness
                    bestFitness = fit_val;
                    bestPosition = SalpPositions(i,:);
                end
            end
        end

        Positions = SalpPositions;
        convergenceCurve(l) = bestFitness;
    end

end

%%---------------------- Initialization Function ----------------------%%
function Positions = initialization(SearchAgents_no,dim,ub,lb)
    Positions = zeros(SearchAgents_no,dim);
    if numel(ub) == 1
        Positions = rand(SearchAgents_no,dim).*(ub-lb) + lb;
    else
        for i = 1:dim
            Positions(:,i) = rand(SearchAgents_no,1) .* (ub(i)-lb(i)) + lb(i);
        end
    end
end
