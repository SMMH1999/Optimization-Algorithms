function [bestFitness, bestPosition, convergenceCurve] = MSO(lb, ub, dim, PopSize, MaxIter, objFun)
    %_______________________________________________________________________________________%
    %  Mirage Search Optimization (MSO)
    %  Developed in MATLAB R2022a
    %
    %  Programmer: Shijie Zhao and Jiahao He
    %  E-mail: zhaoshijie@lntu.edu.cn, lntuhjh@126.com
    %
    %  Description:
    %    MSO is a metaheuristic optimization algorithm inspired by mirage phenomena.
    %    It balances exploration (inferior mirage search) and exploitation (superior mirage search)
    %    to solve global optimization problems.
    %
    %  Inputs:
    %    objFun   - function handle, objective function to minimize
    %    dim      - integer, number of decision variables
    %    lb       - scalar or 1xdim vector, lower bound of variables
    %    ub       - scalar or 1xdim vector, upper bound of variables
    %    PopSize  - integer, population size
    %    MaxIter  - integer, maximum number of iterations
    %
    %  Outputs:
    %    bestFitness      - scalar, best objective function value found
    %    bestPosition     - 1xdim vector, decision variables of the best solution
    %    convergenceCurve - 1xMaxIter vector, best fitness at each iteration
    %
    %  Tunable Parameters:
    %    None explicitly; the algorithm uses adaptive search steps based on iteration
    %_______________________________________________________________________________________%

    %% Initialization
    Positions = initializePopulation(PopSize, dim, ub, lb);
    fitnesses = zeros(PopSize, 1);
    bestFitness = inf;
    bestPosition = zeros(1, dim);
    nfes = 0; % number of function evaluations

    for i = 1:PopSize
        fitnesses(i) = objFun(Positions(i, :));
        nfes = nfes + 1;
        if fitnesses(i) < bestFitness
            bestFitness = fitnesses(i);
            bestPosition = Positions(i, :);
        end
    end

    convergenceCurve = zeros(1, MaxIter);
    iter = 0;

    %% Main Loop
    while iter < MaxIter
        iter = iter + 1;

        %% Superior Mirage Search
        ac = randperm(PopSize-1) + 1;
        cv = ceil((PopSize*(2/3)) * ((MaxIter - nfes + 1) / MaxIter)); % Selection count

        for j = ac(1:cv)
            cosx = zeros(1, dim);
            for k = 1:dim
                h = (bestPosition(k) - Positions(j, k)) * rand();
                cmax = 1;
                h = min(max(h, cmax), 5 * atanh(-(nfes/MaxIter) + 1) + cmax);

                zf = randi(2) * 2 - 3;
                a = rand() * 20;
                b = rand() * (45 - a/2);
                z = randi(2);

                if z == 1
                    C = b + 90;
                    D = 180 - C - a;
                    B = 180 - 2*D;
                    A = 180 - B + a - 90;
                    dx = (sind(B)*h*sind(C)) / (sind(D)*sind(A));
                elseif z == 2 && a < b
                    C = 90 - b;
                    D = 90 + a - b;
                    B = 180 - 2*D;
                    A = 180 - B - a - 90;
                    dx = (sind(B)*h*sind(C)) / (sind(D)*sind(A));
                elseif z == 2 && a > b
                    C = 90 - b;
                    D = 90 - C - a;
                    B = 180 - 2*D;
                    A = 180 - B - 90 + a;
                    dx = (sind(B)*h*sind(C)) / (sind(D)*sind(A));
                else
                    dx = 0;
                end
                cosx(k) = Positions(j, k) + dx * zf;
            end

            % Boundary control
            cosx = min(max(cosx, lb), ub);
            nfes = nfes + 1;
            cosFit = objFun(cosx);

            if cosFit < bestFitness
                bestFitness = cosFit;
                bestPosition = cosx;
            end

            Positions = [Positions; cosx];
            fitnesses = [fitnesses; cosFit];
        end

        % Population renewal
        [fitnesses, idx] = sort(fitnesses);
        Positions = Positions(idx(1:PopSize), :);
        fitnesses = fitnesses(1:PopSize);

        %% Inferior Mirage Search
        for j = 1:PopSize
            hh = bestPosition - Positions(j, :);
            if all(hh == 0)
                hh = 0.05*(randi(2, 1, dim)*2 - 3);
            end

            zf = sign(hh);
            hh = abs(hh .* rand(1, dim));
            gama = rand(1, dim).*90.*((MaxIter - nfes*0.99)/MaxIter);
            amax = atand(1 ./ (2*tand(gama)));
            amin = atand((sind(gama).*cosd(gama)) ./ (1 + sind(gama).^2));
            fai  = (amax - amin).*rand(1, dim) + amin;
            omg  = asind(rand(1, dim).*sind(fai + gama));
            x    = (hh ./ tand(gama)) - ((((hh ./ sind(gama)) - (hh .* sind(fai)) ./ (cosd(fai+gama))).*cosd(omg)) ./ cosd(omg - gama));
            cosx = Positions(j, :) + x.*zf;

            % Boundary control
            cosx = min(max(cosx, lb), ub);
            nfes = nfes + 1;
            cosFit = objFun(cosx);

            if cosFit < bestFitness
                bestFitness = cosFit;
                bestPosition = cosx;
            end

            Positions = [Positions; cosx];
            fitnesses = [fitnesses; cosFit];
        end

        % Population renewal
        [fitnesses, idx] = sort(fitnesses);
        Positions = Positions(idx(1:PopSize), :);
        fitnesses = fitnesses(1:PopSize);

        % Store convergence
        convergenceCurve(iter) = bestFitness;
    end

    %% Local Function: Population Initialization
    function Pop = initializePopulation(nPop, dim, ub, lb)
        if isscalar(ub) && isscalar(lb)
            Pop = rand(nPop, dim) .* (ub - lb) + lb;
        else
            Pop = zeros(nPop, dim);
            for d = 1:dim
                Pop(:, d) = rand(nPop, 1) .* (ub(d) - lb(d)) + lb(d);
            end
        end
    end

end
