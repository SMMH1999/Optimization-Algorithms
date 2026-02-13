function [Best_score, Best_NC, Convergence_curve, t] = NOA(SearchAgents_no, Max_iter, ub, lb, dim, fobj, fhd, Bm)
    %_________________________________________________________________________%
    %  Nutcracker Optimization Algorithm (NOA) - Unified Function              %
    %                                                                          %
    %  Developed in MATLAB R2019A                                              %
    %                                                                          %
    %  Author: Reda Mohamed & Mohamed Abdel-Basset                              %
    %  Reference: Abdel-Basset, M., Mohamed, R., Nutcracker optimizer,         %
    %             Knowledge-Based Systems, in press, 2022,                     %
    %             DOI: https://doi.org/10.1016/j.knosys.2022.110248           %
    %_________________________________________________________________________%
    %
    % Description:
    %   Nutcracker Optimization Algorithm (NOA) is a metaheuristic optimization
    %   method inspired by the foraging and caching behaviors of nutcrackers.
    %
    % Inputs:
    %   SearchAgents_no - integer, number of search agents (nutcrackers)
    %   Max_iter        - integer, maximum number of iterations
    %   ub              - vector or scalar, upper bound(s) of decision variables
    %   lb              - vector or scalar, lower bound(s) of decision variables
    %   dim             - integer, dimension of the problem
    %   fobj            - function handle, objective function to minimize
    %   fhd             - function handle for benchmark (optional, for CEC)
    %   Bm              - integer flag: 1 for benchmark functions, 0 for normal
    %
    % Outputs:
    %   Best_score      - best fitness value found
    %   Best_NC         - best position (solution) found
    %   Convergence_curve - vector storing the best fitness at each iteration
    %   t               - total number of function evaluations
    %
    % Tunable Parameters:
    %   Alpha = 0.05  -> attempts at avoiding local optima
    %   Pa2   = 0.2   -> probability of switching between cache-search and recovery
    %   Prb   = 0.2   -> percentage of exploration in reference points
    %_________________________________________________________________________%


    %%-------------------Parameters--------------------------%%
    Alpha = 0.05;
    Pa2   = 0.2;
    Prb   = 0.2;

    %%-------------------Initialization----------------------%%
    Positions = initializeAgents(SearchAgents_no, dim, ub, lb);
    Lbest     = Positions;
    Best_NC   = zeros(1, dim);
    Convergence_curve = zeros(1, Max_iter);
    LFit = inf(SearchAgents_no,1);
    RP = zeros(2, dim);
    t = 0;

    %%-------------------Initial Fitness Evaluation-----------%%
    NC_Fit = zeros(SearchAgents_no,1);
    for i = 1:SearchAgents_no
        NC_Fit(i) = evaluateFitness(Positions(i,:), fobj, fhd, Bm);
        LFit(i) = NC_Fit(i);
        if NC_Fit(i) < Best_score
            Best_score = NC_Fit(i);
            Best_NC = Positions(i,:);
        end
    end

    %%-------------------Main Optimization Loop----------------%%
    while t < Max_iter
        RL = 0.05 * levy(SearchAgents_no, dim, 1.5);
        l = rand * (1 - t/Max_iter);

        if rand < rand
            a = (t/Max_iter)^(2*1/t);
        else
            a = (1 - t/Max_iter)^(2*(t/Max_iter));
        end

        for i = 1:SearchAgents_no
            if rand < rand
                %% Foraging and Storage Strategy
                mu = selectMu(RL, t, Max_iter);
                cv = randi(SearchAgents_no);
                cv1 = randi(SearchAgents_no);
                Pa1 = (Max_iter - t)/Max_iter;

                if rand < Pa1
                    cv2 = randi(SearchAgents_no);
                    r2 = rand;
                    Positions(i,:) = explorationPhase1(Positions(i,:), Positions(cv,:), Positions(cv1,:), Positions(cv2,:), RL(i,:), mu, l, t, Max_iter, ub, lb, Alpha);
                else
                    Positions(i,:) = exploitationPhase1(Positions(i,:), Positions(cv,:), Positions(cv1,:), Best_NC, RL(i,:), mu);
                end

                %% Enforce bounds
                Positions(i,:) = enforceBounds(Positions(i,:), ub, lb);

                %% Evaluate fitness
                NC_Fit(i) = evaluateFitness(Positions(i,:), fobj, fhd, Bm);

                %% Update local best
                if NC_Fit(i) < LFit(i)
                    LFit(i) = NC_Fit(i);
                    Lbest(i,:) = Positions(i,:);
                else
                    NC_Fit(i) = LFit(i);
                    Positions(i,:) = Lbest(i,:);
                end

                %% Update global best
                if NC_Fit(i) < Best_score
                    Best_score = NC_Fit(i);
                    Best_NC = Positions(i,:);
                end

                t = t + 1;
                if t > Max_iter
                    break;
                end
                Convergence_curve(t) = Best_score;

            else
                %% Cache-Search and Recovery Strategy
                RP = computeReferencePoints(Positions(i,:), Positions, ub, lb, a, Prb);
                if rand < Pa2
                    Positions(i,:) = recoveryStage(Positions(i,:), Best_NC, RP, Positions, ub, lb);
                else
                    Positions(i,:) = cacheSearchStage(Positions(i,:), RP, fobj, fhd, Bm, LFit, Lbest, Best_NC, ub, lb);
                end

                NC_Fit(i) = evaluateFitness(Positions(i,:), fobj, fhd, Bm);

                if NC_Fit(i) < LFit(i)
                    LFit(i) = NC_Fit(i);
                    Lbest(i,:) = Positions(i,:);
                else
                    NC_Fit(i) = LFit(i);
                    Positions(i,:) = Lbest(i,:);
                end

                if NC_Fit(i) < Best_score
                    Best_score = NC_Fit(i);
                    Best_NC = Positions(i,:);
                end

                t = t + 1;
                if t > Max_iter
                    break;
                end
                Convergence_curve(t) = Best_score;
            end
        end
    end

end

%%-------------------Helper Functions--------------------%%

function Positions = initializeAgents(SearchAgents_no, dim, ub, lb)
    if isscalar(ub)
        Positions = rand(SearchAgents_no, dim) * (ub - lb) + lb;
    else
        Positions = rand(SearchAgents_no, dim);
        for i = 1:dim
            Positions(:,i) = Positions(:,i) * (ub(i) - lb(i)) + lb(i);
        end
    end
end

function fit = evaluateFitness(pos, fobj, fhd, Bm)
    if Bm >= 1
        fit = feval(fhd, pos', fobj);
    else
        fit = fobj(pos);
    end
end

function z = levy(n, m, beta)
    num = gamma(1 + beta) * sin(pi*beta/2);
    den = gamma((1 + beta)/2) * beta * 2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta);
    u = randn(n,m) * sigma_u;
    v = randn(n,m);
    z = u ./ abs(v).^(1/beta);
end

function mu = selectMu(RL, t, Max_iter)
    r = rand;
    if r < rand
        mu = rand;
    elseif r < rand
        mu = randn;
    else
        mu = RL(1,1);
    end
end

function pos = explorationPhase1(pos, cv_pos, cv1_pos, cv2_pos, RLi, mu, l, t, Max_iter, ub, lb, Alpha)
    dim = length(pos);
    r2 = rand;
    for j = 1:dim
        if t < Max_iter/2
            if rand > rand
                pos(j) = (mean([cv_pos(j), cv1_pos(j)])) + RLi(j) * (cv_pos(j) - cv1_pos(j)) + mu * (r2*r2*ub(j) - lb(j));
            end
        else
            if rand > rand
                pos(j) = cv2_pos(j) + mu * (cv_pos(j) - cv1_pos(j)) + mu * (rand < Alpha) * (r2*r2*ub(j) - lb(j));
            end
        end
    end
end

function pos = exploitationPhase1(pos, cv_pos, cv1_pos, Best_NC, RLi, mu)
    dim = length(pos);
    for j = 1:dim
        r1 = rand;
        if rand < rand
            pos(j) = pos(j) + mu * abs(RLi(j)) * (Best_NC(j) - pos(j)) + r1 * (cv_pos(j) - cv1_pos(j));
        elseif rand < rand
            if rand > rand
                pos(j) = Best_NC(j) + mu * (cv_pos(j) - cv1_pos(j));
            end
        else
            pos(j) = Best_NC(j) * abs(rand);
        end
    end
end

function pos = enforceBounds(pos, ub, lb)
    pos = min(max(pos, lb), ub);
end

function RP = computeReferencePoints(pos, Positions, ub, lb, a, Prb)
    dim = length(pos);
    RP = zeros(2, dim);
    cv = randi(size(Positions,1));
    cv1 = randi(size(Positions,1));
    ang = pi*rand;
    for j = 1:dim
        if ang ~= pi/2
            RP(1,j) = pos(j) + a*cos(ang)*(Positions(cv,j) - Positions(cv1,j));
            RP(2,j) = pos(j) + a*cos(ang)*((ub(j)-lb(j)) + lb(j))*(rand < Prb);
        else
            RP(1,j) = pos(j) + a*cos(ang)*(Positions(cv,j) - Positions(cv1,j)) + a*RP(randi(2), j);
            RP(2,j) = pos(j) + a*cos(ang)*((ub(j)-lb(j))*rand + lb(j)) + a*RP(randi(2), j)*(rand < Prb);
        end
    end
    RP = min(max(RP, lb), ub);
end

function pos = recoveryStage(pos, Best_NC, RP, Positions, ub, lb)
    cv = randi(size(Positions,1));
    for j = 1:length(pos)
        if rand > rand
            pos(j) = pos(j) + rand*(Best_NC(j) - pos(j)) + rand*(RP(1,j) - Positions(cv,j));
        else
            pos(j) = pos(j) + rand*(Best_NC(j) - pos(j)) + rand*(RP(2,j) - Positions(cv,j));
        end
    end
    pos = min(max(pos, lb), ub);
end

function pos = cacheSearchStage(pos, RP, fobj, fhd, Bm, LFit, Lbest, Best_NC, ub, lb)
    NC_Fit1 = evaluateFitness(RP(1,:), fobj, fhd, Bm);
    NC_Fit2 = evaluateFitness(RP(2,:), fobj, fhd, Bm);
    NC_Fit = evaluateFitness(pos, fobj, fhd, Bm);

    if NC_Fit2 < NC_Fit1 && NC_Fit2 < NC_Fit
        pos = RP(2,:);
    elseif NC_Fit1 < NC_Fit2 && NC_Fit1 < NC_Fit
        pos = RP(1,:);
    end

    pos = min(max(pos, lb), ub);
end
