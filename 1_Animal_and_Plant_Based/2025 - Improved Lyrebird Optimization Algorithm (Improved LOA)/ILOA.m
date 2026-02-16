function [bestFitness, bestPosition, convergenceCurve] = ILOA(lb, ub, dim, nPop, maxItr, objFun)
%%----------------------------------------------------------------------
% Improved Lyrebird Optimization Algorithm (ILOA)
% Author: A. Khaled et al., inspired by Weight Minimization problem demo
%
% Inputs:
%   lb        : Lower bounds (1 x dim or scalar)
%   ub        : Upper bounds (1 x dim or scalar)
%   dim       : Number of dimensions
%   nPop      : Population size (number of lyrebirds)
%   maxItr    : Maximum number of iterations
%   objFun    : Handle to the objective function
%
% Outputs:
%   bestFitness       : Best fitness value found
%   bestPosition      : Position vector of best solution
%   convergenceCurve  : Best fitness at each iteration
%
% Tunable Parameters:
%   P_m     : Probability of memory-based search (mimicry)
%   gamma   : Repulsion coefficient for territory protection
%----------------------------------------------------------------------

% Parameters
P_m = 0.2;      % Mimicry probability
gamma = 0.1;    % Territory protection coefficient

% Ensure bounds are row vectors
if isscalar(lb), lb = lb*ones(1,dim); end
if isscalar(ub), ub = ub*ones(1,dim); end

% Initialization
X = lb + rand(nPop, dim) .* (ub - lb);      % Population
fitness = zeros(nPop,1);
for i = 1:nPop
    fitness(i) = objFun(X(i,:));
end
memory = X;                                 % Memory for mimicry

% Best initialization
[bestFitness, bestIdx] = min(fitness);
bestPosition = X(bestIdx,:);

% Convergence tracking
convergenceCurve = zeros(1, maxItr);

% Main loop
for t = 1:maxItr
    D_max = max(pdist(X)); % Maximum observed distance
    
    for i = 1:nPop
        rp = rand(); % Decide phase
        
        if rp <= 0.5 % Phase 1: Escaping (Exploration)
            % Identify safe areas
            safeAreas = X(fitness < fitness(i), :);
            if ~isempty(safeAreas)
                safeArea = safeAreas(randi(size(safeAreas,1)), :);
                r1 = rand(1,dim);
                I = randi([1,2],1,dim);
                X_new = X(i,:) + r1 .* (safeArea - I .* X(i,:));
            else
                X_new = X(i,:);
            end
            
            % Mimicry mechanism (memory-based)
            if rand() <= P_m*t/maxItr
                M_r1 = memory(randi(nPop), :);
                X_new = bestPosition + rand(1,dim) .* (M_r1 - X(i,:));
            end
            
        else % Phase 2: Hiding (Exploitation)
            r2 = rand(1,dim);
            X_new = X(i,:) + (1 - 2*r2) .* ((ub - lb) / t);
            
            % Territory protection (diversity)
            D_i = mean(pdist2(X(i,:), X));
            if D_i < gamma*D_max
                X_j = X(randi(nPop), :);
                X_new = bestPosition + gamma*rand(1,dim) .* (X(i,:) - X_j);
            end
        end
        
        % Enforce boundaries
        X_new = max(min(X_new, ub), lb);
        
        % Evaluate new fitness
        newFitness = objFun(X_new);
        
        % Update if improved
        if newFitness < fitness(i)
            X(i,:) = X_new;
            fitness(i) = newFitness;
            memory(i,:) = X_new; % Update memory
        end
    end
    
    % Update global best
    [currentBest, bestIdx] = min(fitness);
    if currentBest < bestFitness
        bestFitness = currentBest;
        bestPosition = X(bestIdx,:);
    end
    
    % Record convergence
    convergenceCurve(t) = bestFitness;
end

end
