function [bestFitness, bestPosition, convergenceCurve] = FATA(lb, ub, dim, nPop, maxItr, objFun)
    %=================================================================================
    % FATA: Fata Morgana Optimization Algorithm
    % Version: 1.0 (Refactored for benchmark framework)
    %
    % Author: Ailiang Qi, Dong Zhao, Ali Asghar Heidari, Lei Liu, Yi Chen, Huiling Chen
    % Reference: FATA: An Efficient Optimization Method Based on Geophysics, Neurocomputing, 2024
    %
    % Description:
    %   FATA is a nature-inspired optimization algorithm simulating light refraction
    %   and mirage effects to explore the search space effectively.
    %
    % Inputs:
    %   lb       - lower bounds (scalar or 1 x dim vector)
    %   ub       - upper bounds (scalar or 1 x dim vector)
    %   dim      - number of decision variables
    %   nPop     - population size (number of agents)
    %   maxItr   - maximum number of iterations
    %   objFun   - handle of objective function to minimize
    %
    % Outputs:
    %   bestFitness       - best objective function value found
    %   bestPosition      - decision variable vector corresponding to bestFitness
    %   convergenceCurve  - vector of bestFitness at each iteration
    %
    % Tunable Parameters:
    %   arf - reflection factor for light internal reflection (default 0.2)
    %=================================================================================

    %------------------------ Initialization ------------------------%
    arf = 0.2;                      % Light total internal reflection factor
    bestFitness = inf;               % Initialize best fitness
    bestPosition = zeros(1, dim);   % Initialize best position
    convergenceCurve = zeros(1, maxItr);

    % Ensure bounds are vectors
    lb = ones(1, dim) * lb;
    ub = ones(1, dim) * ub;

    % Initialize population
    Population = initialization(nPop, dim, ub, lb);
    fitness = inf(nPop, 1);

    % Initialize population quality integrals
    worstInte = 0;
    bestInte = inf;

    %------------------------ Main Loop ------------------------%
    for it = 1:maxItr
        % Boundary check & evaluate fitness
        for i = 1:nPop
            Population(i,:) = min(max(Population(i,:), lb), ub);
            fitness(i) = objFun(Population(i,:));
            if fitness(i) < bestFitness
                bestFitness = fitness(i);
                bestPosition = Population(i,:);
            end
        end

        % Sort fitness
        [sortedFit, ~] = sort(fitness);
        worstFitness = sortedFit(end);
        bestFit = sortedFit(1);

        % Mirage light filtering principle (population quality factor)
        Integral = cumtrapz(sortedFit);
        worstInte = max(worstInte, Integral(end));
        bestInte = min(bestInte, Integral(end));
        IP = (Integral(end) - worstInte) / (bestInte - worstInte + eps);

        % Parameters for refraction
        a = tan(-(it/maxItr) + 1);
        b = 1 / tan(-(it/maxItr) + 1);

        % Update positions
        for i = 1:nPop
            Para1 = a*rand(1,dim) - a*rand(1,dim);
            Para2 = b*rand(1,dim) - b*rand(1,dim);
            p = (fitness(i) - worstFitness) / (bestFitness - worstFitness + eps);

            if rand > IP
                Population(i,:) = (ub - lb).*rand(1, dim) + lb;
            else
                for j = 1:dim
                    idx = randi(nPop);
                    if rand < p
                        % Light refraction (first phase)
                        Population(i,j) = bestPosition(j) + Population(i,j) * Para1(j);
                    else
                        % Light refraction (second phase)
                        Population(i,j) = Population(idx,j) + Para2(j) * Population(i,j);
                        % Total internal reflection
                        Population(i,j) = (0.5*(arf+1)*(lb(j)+ub(j)) - arf*Population(i,j));
                    end
                end
            end
        end

        convergenceCurve(it) = bestFitness;
    end

end

%=================== Helper Function ===================%
function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Positions = zeros(SearchAgents_no, dim);
    for i = 1:dim
        Positions(:,i) = rand(SearchAgents_no,1) .* (ub(i) - lb(i)) + lb(i);
    end
end
