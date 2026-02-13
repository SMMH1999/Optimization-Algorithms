function [Best_Score, BestFit, Convergence_curve] = SCSO(lb, ub, dim, SearchAgents_no, Max_iter, fobj)
    %% Sand Cat Optimization Algorithm (SCSO)
    % Author: Amir Seyyedabbasi
    % e-Mail: amir.seyedabbasi@gmail.com
    % Main paper: A. Seyyedabbasi, F. Kiani
    % DOI: https://doi.org/10.1007/s00366-022-01604-x
    %
    % Description:
    %   This function implements the Sand Cat Optimization Algorithm (SCSO)
    %   for continuous optimization problems.
    %
    % Inputs:
    %   SearchAgents_no - (int) Number of search agents (population size)
    %   Max_iter        - (int) Maximum number of iterations
    %   lb              - (scalar or 1xdim vector) Lower bounds of variables
    %   ub              - (scalar or 1xdim vector) Upper bounds of variables
    %   dim             - (int) Number of variables (dimensions)
    %   fobj            - (function handle) Objective function to minimize
    %
    % Outputs:
    %   Best_Score       - (double) Best objective function value found
    %   BestFit          - (1xdim double) Position vector of the best solution
    %   Convergence_curve - (1xMax_iter double) Best score at each iteration

    %% Initialization
    BestFit = zeros(1, dim);           % Best position
    Best_Score = inf;                  % Best score (to minimize)
    Positions = initializePopulation(SearchAgents_no, dim, ub, lb);  % Initial population
    Convergence_curve = zeros(1, Max_iter);
    t = 0;                             % Iteration counter
    p = 1:360;                         % Roulette wheel selection probabilities

    %% Main Loop
    while t < Max_iter
        %% Fitness Evaluation and Best Update
        for i = 1:size(Positions, 1)
            % Boundary control
            Flag4ub = Positions(i,:) > ub;
            Flag4lb = Positions(i,:) < lb;
            Positions(i,:) = (Positions(i,:) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;

            % Fitness calculation
            fitness = fobj(Positions(i,:));

            % Update best solution
            if fitness < Best_Score
                Best_Score = fitness;
                BestFit = Positions(i,:);
            end
        end

        %% Position Update
        S = 2;                                 % Maximum sensitivity range
        rg = S - (S * t / Max_iter);           % Guides R

        for i = 1:size(Positions,1)
            r = rand * rg;
            R = ((2 * rg) * rand) - rg;        % Transition control

            for j = 1:size(Positions,2)
                teta = RouletteWheelSelection(p);

                if (-1 <= R) && (R <= 1)
                    Rand_position = abs(rand * BestFit(j) - Positions(i,j));
                    Positions(i,j) = BestFit(j) - r * Rand_position * cos(teta);
                else
                    cp = floor(SearchAgents_no * rand() + 1);
                    CandidatePosition = Positions(cp,:);
                    Positions(i,j) = r * (CandidatePosition(j) - rand * Positions(i,j));
                end
            end
        end

        %% Convergence curve update
        t = t + 1;
        Convergence_curve(t) = Best_Score;
    end

end

%% ---------------- Helper Functions ---------------- %%
function j = RouletteWheelSelection(P)
    r = rand;
    P = P / sum(P);
    C = cumsum(P);
    j = find(r <= C, 1, 'first');
end

function Positions = initializePopulation(SearchAgents_no, dim, ub, lb)
    Boundary_no = numel(ub); % Number of boundaries

    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) * (ub - lb) + lb;
    else
        Positions = zeros(SearchAgents_no, dim);
        for i = 1:dim
            Positions(:,i) = rand(SearchAgents_no,1) * (ub(i) - lb(i)) + lb(i);
        end
    end
end
