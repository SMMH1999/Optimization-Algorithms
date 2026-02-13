function [Top_predator_fit,Top_predator_pos,Convergence_curve]=MPA(lb,ub,dim,SearchAgents_no,Max_iter,fobj)
    %--------------------------------------------------------------------------
    % Marine Predators Algorithm (MPA)
    %--------------------------------------------------------------------------
    % Author: Afshin Faramarzi & Seyedali Mirjalili
    % Paper: A. Faramarzi et al., Marine Predators Algorithm: A Nature-inspired Metaheuristic,
    %        Expert Systems with Applications, 2020
    %--------------------------------------------------------------------------
    % INPUTS:
    %   SearchAgents_no : int    : Number of search agents (population size)
    %   Max_iter        : int    : Maximum number of iterations
    %   lb              : float  : Lower bound (scalar or 1xdim vector)
    %   ub              : float  : Upper bound (scalar or 1xdim vector)
    %   dim             : int    : Number of decision variables
    %   fobj            : func   : Objective function handle
    %
    % OUTPUTS:
    %   Top_predator_fit : float : Best objective value found
    %   Top_predator_pos : 1xdim : Position of best solution
    %   Convergence_curve: 1xMax_iter : Best fitness value at each iteration
    %
    % Tunable parameters:
    %   FADs : float : Probability of Fish Aggregating Devices effect (default 0.2)
    %   P    : float : Predation parameter (default 0.5)
    %--------------------------------------------------------------------------
    % Initialization
    Top_predator_pos = zeros(1,dim);
    Top_predator_fit = inf;
    Convergence_curve = zeros(1,Max_iter);
    stepsize = zeros(SearchAgents_no,dim);
    fitness = inf(SearchAgents_no,1);

    Prey = initialization(SearchAgents_no,dim,ub,lb);
    Xmin = repmat(lb, SearchAgents_no,1);
    Xmax = repmat(ub, SearchAgents_no,1);

    Iter = 0;
    FADs = 0.2;
    P = 0.5;

    %---------------------- Main Loop ----------------------%
    while Iter < Max_iter
        %------------------- Boundary check & Fitness evaluation -----------------
        for i=1:SearchAgents_no
            Prey(i,:) = max(min(Prey(i,:), ub), lb);
            fitness(i) = fobj(Prey(i,:));
            if fitness(i) < Top_predator_fit
                Top_predator_fit = fitness(i);
                Top_predator_pos = Prey(i,:);
            end
        end

        %------------------- Marine Memory saving -------------------
        if Iter == 0
            fit_old = fitness;
            Prey_old = Prey;
        end
        Inx = (fit_old < fitness);
        Indx = repmat(Inx,1,dim);
        Prey = Indx.*Prey_old + ~Indx.*Prey;
        fitness = Inx.*fit_old + ~Inx.*fitness;
        fit_old = fitness;
        Prey_old = Prey;

        %------------------- Elite & Random Matrices -------------------
        Elite = repmat(Top_predator_pos, SearchAgents_no,1);
        CF = (1 - Iter/Max_iter)^(2*Iter/Max_iter);
        RL = 0.05*levy(SearchAgents_no,dim,1.5);
        RB = randn(SearchAgents_no,dim);

        %------------------- Update Prey positions -------------------
        for i=1:SearchAgents_no
            for j=1:dim
                R = rand();
                if Iter < Max_iter/3
                    % Phase 1
                    stepsize(i,j) = RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));
                    Prey(i,j) = Prey(i,j) + P*R*stepsize(i,j);
                elseif Iter < 2*Max_iter/3
                    % Phase 2
                    if i > SearchAgents_no/2
                        stepsize(i,j) = RB(i,j)*(RB(i,j)*Elite(i,j)-Prey(i,j));
                        Prey(i,j) = Elite(i,j) + P*CF*stepsize(i,j);
                    else
                        stepsize(i,j) = RL(i,j)*(Elite(i,j)-RL(i,j)*Prey(i,j));
                        Prey(i,j) = Prey(i,j) + P*R*stepsize(i,j);
                    end
                else
                    % Phase 3
                    stepsize(i,j) = RL(i,j)*(RL(i,j)*Elite(i,j)-Prey(i,j));
                    Prey(i,j) = Elite(i,j) + P*CF*stepsize(i,j);
                end
            end
        end

        %------------------- Boundary check & Fitness evaluation -----------------
        for i=1:SearchAgents_no
            Prey(i,:) = max(min(Prey(i,:), ub), lb);
            fitness(i) = fobj(Prey(i,:));
            if fitness(i) < Top_predator_fit
                Top_predator_fit = fitness(i);
                Top_predator_pos = Prey(i,:);
            end
        end

        %------------------- Marine Memory saving -------------------
        Inx = (fit_old < fitness);
        Indx = repmat(Inx,1,dim);
        Prey = Indx.*Prey_old + ~Indx.*Prey;
        fitness = Inx.*fit_old + ~Inx.*fitness;
        fit_old = fitness;
        Prey_old = Prey;

        %---------- Eddy formation and FADs’ effect ----------------
        if rand() < FADs
            U = rand(SearchAgents_no,dim) < FADs;
            Prey = Prey + CF*((Xmin + rand(SearchAgents_no,dim).*(Xmax-Xmin)).*U);
        else
            r = rand();  Rs = SearchAgents_no;
            stepsize = (FADs*(1-r) + r)*(Prey(randperm(Rs),:) - Prey(randperm(Rs),:));
            Prey = Prey + stepsize;
        end

        Iter = Iter + 1;
        Convergence_curve(Iter) = Top_predator_fit;
    end
end

%---------------------- Levy Flight ----------------------%
function z = levy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta);
    u = random('Normal',0,sigma_u,n,m);
    v = random('Normal',0,1,n,m);
    z = u ./ (abs(v).^(1/beta));
end

%---------------------- Initialization ----------------------%
function Positions = initialization(SearchAgents_no,dim,ub,lb)
    if numel(ub)==1
        Positions = rand(SearchAgents_no,dim).*(ub-lb) + lb;
    else
        Positions = zeros(SearchAgents_no,dim);
        for i=1:dim
            Positions(:,i) = rand(SearchAgents_no,1).*(ub(i)-lb(i)) + lb(i);
        end
    end
end
