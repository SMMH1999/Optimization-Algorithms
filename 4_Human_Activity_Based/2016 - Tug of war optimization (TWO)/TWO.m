function [bestFitness, bestPosition, convergenceCurve] = TWO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Tug of War Optimization (TWO) - Pairwise Competition Version
    % =========================================================================
    % Description:
    % Teams are sorted by fitness. Weights are assigned inversely proportional
    % to fitness. Weaker teams move toward stronger teams using pairwise
    % displacement mechanism.
    %
    % Inputs:
    % lb, ub   : Lower/Upper bounds (scalar or vector)
    % dim      : Dimension
    % nPop     : Population size
    % maxItr   : Maximum iterations
    % objFun   : Objective function handle (minimization)
    %
    % Outputs:
    % bestFitness      : Best objective value
    % bestPosition     : Best solution
    % convergenceCurve : Best fitness at each iteration
    % =========================================================================

    %% Bound handling
    if numel(lb)==1, lb=lb*ones(1,dim); end
    if numel(ub)==1, ub=ub*ones(1,dim); end

    %% Initialization
    solutions=rand(nPop,dim).*(ub-lb)+lb;
    fitness=zeros(nPop,1);

    for i=1:nPop
        fitness(i)=objFun(solutions(i,:));
    end

    [bestFitness,idx]=min(fitness);
    bestPosition=solutions(idx,:);
    convergenceCurve=zeros(maxItr,1);

    %% Main loop
    for t=1:maxItr

        [fitness,sortedIdx]=sort(fitness);
        solutions=solutions(sortedIdx,:);

        weights=1./(1+fitness);

        for i=1:nPop
            for j=1:nPop
                if i~=j && weights(i)<weights(j)
                    displacement=rand*(solutions(j,:)-solutions(i,:));
                    solutions(i,:)=solutions(i,:)+displacement;
                end
            end

            solutions(i,:)=max(solutions(i,:),lb);
            solutions(i,:)=min(solutions(i,:),ub);
            fitness(i)=objFun(solutions(i,:));
        end

        [currentBest,idx]=min(fitness);
        if currentBest<bestFitness
            bestFitness=currentBest;
            bestPosition=solutions(idx,:);
        end

        convergenceCurve(t)=bestFitness;
    end
end
