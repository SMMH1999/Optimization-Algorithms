function [bestFitness, bestPosition, convergenceCurve] = LA(lb, ub, dim, nPop, maxItr, objFun)
    %--------------------------------------------------------------------------
    % Lichtenberg Algorithm (LA) - Core Optimization Only
    %
    % Author: Original: Qamar Askari
    %
    % Description:
    %   Population-based metaheuristic optimization algorithm inspired by
    %   Lichtenberg figures. This version contains only the core algorithm,
    %   no visualization.
    %
    % Inputs:
    %   lb        : Lower bound (scalar or 1xd vector)
    %   ub        : Upper bound (scalar or 1xd vector)
    %   dim       : Number of decision variables
    %   nPop      : Population size
    %   maxItr    : Maximum number of iterations
    %   objFun    : Function handle for objective function (supports constraints)
    %
    % Outputs:
    %   bestFitness       : Best objective function value found
    %   bestPosition      : Best solution vector found
    %   convergenceCurve  : Best fitness value at each iteration
    %--------------------------------------------------------------------------
    % Parameters for Lichtenberg Figure (tunable)
    Np = 50;    % Number of points in figure
    Rc = 5;     % Radius constant
    S_c = 0.5;  % Sticking probability
    M = 0;      % Figure mode: 0=random, 1,2=dynamic
    ref = 0;    % Reference scaling (0=off)

    % Ensure lb and ub are vectors
    if isscalar(lb), lb = repmat(lb,1,dim); end
    if isscalar(ub), ub = repmat(ub,1,dim); end

    %--------------------------------------------------------------------------
    % Initialize Population
    Individuals = lb + (ub-lb).*rand(nPop,dim);
    Fitness = zeros(nPop,1);
    for i=1:nPop
        Fitness(i) = objFun(Individuals(i,:));
    end

    % Best initialization
    [bestFitness,I] = min(Fitness);
    bestPosition = Individuals(I,:);
    convergenceCurve = zeros(maxItr,1);

    % Initial figure points
    S = Individuals;
    if M==0 || M==1
        LF = LA_figure(dim,Np,Rc,S_c,M);
    end

    %--------------------------------------------------------------------------
    % Main Optimization Loop
    for t=1:maxItr
        x_start = bestPosition;
        if M==2
            LF = LA_figure(dim,Np,Rc,S_c,M);
        end

        scale_factor = 1.2*rand;
        X = LA_points(LF,lb,ub,x_start,scale_factor,dim);

        if ref ~= 0
            X_local = LA_points(LF,lb*ref,ub*ref,x_start,scale_factor,dim);
        end

        for i=1:nPop
            if ref ~= 0
                pop1 = round(0.4*nPop);
                pop2 = nPop - pop1;
                S_global = X(randperm(size(X,1),pop2),:);
                S_ref = X_local(randperm(size(X_local,1),pop1),:);
                S = [S_global; S_ref];
            else
                S = X(randperm(size(X,1),nPop),:);
            end

            % Boundary check
            S(i,:) = bound_check(S(i,:), lb, ub);

            % Evaluate fitness
            Fnew = objFun(S(i,:));

            % Update individual
            if Fnew <= Fitness(i)
                Individuals(i,:) = S(i,:);
                Fitness(i) = Fnew;
            end

            % Update global best
            if Fnew <= bestFitness
                bestPosition = S(i,:);
                bestFitness = Fnew;
            end
        end

        % Store convergence
        convergenceCurve(t) = bestFitness;
    end
end

%--------------------------------------------------------------------------
function s = bound_check(s,LB,UB)
    s = max(min(s,UB),LB);
end

%--------------------------------------------------------------------------
function [map] = LA_figure(d,Np,Rc,S,M)
    % Core Lichtenberg Figure generation (simplified placeholder)
    if M==0
        map = zeros(2*Rc+4,2*Rc+4, max(1,d));
    else
        map = rand(Np,d)*Rc; % Placeholder: points generation
    end
end

%--------------------------------------------------------------------------
function X = LA_points(K,LB,UB,x0,scale_factor,d)
    % Core transformation of figure points into new population
    X = K;
    for i=1:d
        scale = scale_factor*(UB(i)-LB(i))/(max(X(:,i))-min(X(:,i))+eps);
        X(:,i) = scale*X(:,i);
        delta = mean(X(:,i)) - x0(i);
        X(:,i) = X(:,i) - delta;
    end
end
