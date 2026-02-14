function [bestFitness, bestPosition, convergenceCurve] = HBO(lb, ub, dim, nPop, maxItr, objFun)
    % =========================================================================
    % Hannibal's Battle Optimization (HBO)
    % =========================================================================
    % Author: (Original HBO Authors)
    % Refactored for benchmark format (Instructor/CEC style)
    %
    % Description:
    % Hannibal's Battle Optimization (HBO) is a population-based metaheuristic
    % inspired by Hannibal’s battle strategies. The population is divided into
    % Romans and Carthaginians, which are further separated into strategic
    % groups (left, right, center). The algorithm simulates battle phases:
    % early attack, strategic retreat, and encirclement.
    %
    % =========================================================================
    % Inputs:
    % lb        : Lower bound (scalar or 1×dim vector)
    % ub        : Upper bound (scalar or 1×dim vector)
    % dim       : Problem dimension
    % nPop      : Total population size
    % maxItr    : Maximum number of iterations
    % objFun    : Objective function handle (minimization)
    %
    % Outputs:
    % bestFitness      : Best objective function value found
    % bestPosition     : Best decision variable vector
    % convergenceCurve : Best fitness value at each iteration
    % =========================================================================

    %% --------------------------- Initialization -----------------------------

    if isscalar(lb)
        lb = lb * ones(1, dim);
    end
    if isscalar(ub)
        ub = ub * ones(1, dim);
    end

    bestFitness = inf;
    bestPosition = zeros(1, dim);
    Rom_Gen = zeros(1, dim);
    Rom_score = inf;

    convergenceCurve = zeros(1, maxItr);

    roman_rate = 0.63;
    romans_no = round(nPop * roman_rate);
    cartagos_no = nPop - romans_no;

    romans = initialization(romans_no, dim, ub, lb);
    cartagos = initialization(cartagos_no, dim, ub, lb);

    romans_fitness = zeros(romans_no,1);
    for i = 1:romans_no
        romans_fitness(i) = objFun(romans(i,:));
        if romans_fitness(i) < bestFitness
            bestFitness = romans_fitness(i);
            bestPosition = romans(i,:);
        end
    end

    cartagos_fitness = zeros(cartagos_no,1);
    for i = 1:cartagos_no
        cartagos_fitness(i) = objFun(cartagos(i,:));
        if cartagos_fitness(i) < bestFitness
            bestFitness = cartagos_fitness(i);
            bestPosition = cartagos(i,:);
        end
    end

    %% ---------------------- Divide Strategic Groups -------------------------

    [~, idxR] = sort(romans_fitness);
    romans = romans(idxR,:);
    p3 = max(1, round(0.03 * romans_no));
    left_romans   = romans(1:p3,:);
    right_romans  = romans(p3+1:min(2*p3,romans_no),:);
    center_romans = romans(min(2*p3+1,romans_no):end,:);

    [~, idxC] = sort(cartagos_fitness,'descend');
    cartagos = cartagos(idxC,:);
    p12 = max(1, round(0.12 * cartagos_no));
    left_cartagos   = cartagos(1:p12,:);
    right_cartagos  = cartagos(p12+1:min(2*p12,cartagos_no),:);
    center_cartagos = cartagos(min(2*p12+1,cartagos_no):end,:);

    %% ----------------------------- Main Loop --------------------------------

    for t = 1:maxItr

        if t <= 0.33*maxItr
            phase = 1;
        elseif t <= 0.66*maxItr
            phase = 2;
        else
            phase = 3;
        end

        switch phase
            case 1
                [right_cartagos,left_romans] = right_attack(right_cartagos,left_romans,Rom_Gen,bestPosition);
                [left_cartagos,right_romans] = left_attack(left_cartagos,right_romans,Rom_Gen,bestPosition);
                [center_cartagos,center_romans] = center_attack(center_cartagos,center_romans,Rom_Gen,bestPosition);

            case 2
                [right_cartagos,left_romans] = right_attack(right_cartagos,left_romans,Rom_Gen,bestPosition);
                [left_cartagos,left_romans] = left_attack(left_cartagos,left_romans,Rom_Gen,bestPosition);
                right_romans = right_romans + randn(size(right_romans)).*(ub-lb);
                [center_cartagos,center_romans] = center_attack(center_cartagos,center_romans,Rom_Gen,bestPosition);

            case 3
                [right_cartagos,center_romans] = right_attack(right_cartagos,center_romans,Rom_Gen,bestPosition);
                left_romans = left_romans + randn(size(left_romans)).*(ub-lb);
                [left_cartagos,center_romans] = left_attack(left_cartagos,center_romans,Rom_Gen,bestPosition);
                right_romans = right_romans + randn(size(right_romans)).*(ub-lb);
                [center_cartagos,center_romans] = center_attack(center_cartagos,center_romans,Rom_Gen,bestPosition);
        end

        groups = {left_romans,right_romans,center_romans,...
            left_cartagos,right_cartagos,center_cartagos};

        for g = 1:length(groups)
            pop = groups{g};
            for i = 1:size(pop,1)
                pop(i,:) = enforceBounds(pop(i,:),ub,lb);
                fitness = objFun(pop(i,:));

                if fitness < bestFitness
                    bestFitness = fitness;
                    bestPosition = pop(i,:);
                elseif fitness < Rom_score
                    Rom_score = fitness;
                    Rom_Gen = pop(i,:);
                end
            end
            groups{g} = pop;
        end

        [left_romans,right_romans,center_romans,...
            left_cartagos,right_cartagos,center_cartagos] = deal(groups{:});

        convergenceCurve(t) = bestFitness;
    end

end

%% ============================ Local Functions ===========================

function Positions = initialization(n,dim,ub,lb)
    Positions = zeros(n,dim);
    for i=1:dim
        Positions(:,i)=rand(n,1).*(ub(i)-lb(i))+lb(i);
    end
end

function position = enforceBounds(position, ub, lb)
    position = min(max(position,lb),ub);
end

function [newR,newL] = left_attack(posR,posL,General,Hannibal)
    rvR = randn(size(posR))*0.5;
    rvL = randn(size(posL))*0.5;
    newR = 2*mean(posL) - rvR.*posR;
    newL = 2*mean(posR) - rvL.*posL;
    newR = 0.94*parallaxe(newR,Hannibal);
    newL = 0.94*parallaxe(newL,General);
end

function [newR,newL] = right_attack(posR,posL,General,Hannibal)
    plane = (mean(posR)+mean(posL))/2;
    rvR = randn(size(posR))*0.5;
    rvL = randn(size(posL))*0.5;
    newR = plane + rvR.*(plane-posR);
    newL = plane + rvL.*(plane-posL);
    newR = 0.94*parallaxe(newR,Hannibal);
    newL = 0.94*parallaxe(newL,General);
end

function [newR,newL] = center_attack(posR,posL,General,Hannibal)
    c=0.94;
    rvR = randn(size(posR))*c*rand;
    rvL = randn(size(posL))*c*rand;
    newR = posR + rvR.*(General-posR);
    newL = posL + rvL.*(Hannibal-posL);
end

function modified = parallaxe(pos,Gen)
    N=size(pos,1);
    if N<2
        modified=pos; return;
    end
    if mod(N,2)~=0
        couples=reshape(randperm(N-1),[],2);
    else
        couples=reshape(randperm(N),[],2);
    end
    modified=pos;
    for k=1:size(couples,1)
        c=(pos(couples(k,1),:)+pos(couples(k,2),:))/2;
        d=norm(c-Gen);
        for j=1:2
            idx=couples(k,j);
            dir=(Gen-pos(idx,:));
            if norm(dir)==0, continue; end
            dir=dir/norm(dir);
            step=d*(0.5+rand()*2);
            modified(idx,:)=pos(idx,:)+step*dir;
        end
    end
end
