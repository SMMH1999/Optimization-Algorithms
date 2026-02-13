function [bestFitness, bestPosition, convergenceCurve] = CCO(lb, ub, dim, nPop, maxItr, objFun)
    %_________________________________________________________________________________
    %  Cuckoo Catfish Optimizer (CCO)
    %_________________________________________________________________________________
    %  Original Paper:
    %  Tian-Lei Wang, Shao-Wei Gu, Ren-Ju Liu, Le-Qing Chen, Zhu Wang, Zhi-Qiang Zeng
    %  "Cuckoo Catfish Optimizer: A New Meta-Heuristic Optimization Algorithm"
    %  Artificial Intelligence Review
    %
    %  Refactored to Benchmark / Instructor Format
    %
    %  INPUTS:
    %    lb        : Lower bound (scalar or 1×dim vector)
    %    ub        : Upper bound (scalar or 1×dim vector)
    %    dim       : Problem dimension
    %    nPop      : Population size
    %    maxItr    : Maximum number of iterations
    %    objFun    : Objective function handle (minimization)
    %
    %  OUTPUTS:
    %    bestFitness      : Best objective function value found
    %    bestPosition     : Best solution vector found (1×dim)
    %    convergenceCurve : Best fitness value at each iteration (maxItr×1)
    %
    %  Tunable Parameters (as in original paper):
    %    alpha = 1.34
    %    beta  = 0.3
    %
    %  Notes:
    %    - Original algorithm logic strictly preserved.
    %    - Boundary handling and Lévy flight implemented internally.
    %    - Fully compatible with benchmark optimization frameworks.
    %_________________________________________________________________________________

    %% Bound handling
    if isscalar(lb), lb = lb * ones(1, dim); end
    if isscalar(ub), ub = ub * ones(1, dim); end

    alpha = 1.34;
    beta  = 0.3;

    %% Initialization
    PopPos = rand(nPop, dim) .* (ub - lb) + lb;
    PopFit = zeros(nPop,1);

    bestFitness  = inf;
    bestPosition = zeros(1,dim);

    x = zeros(nPop,1);
    y = zeros(nPop,1);

    for i = 1:nPop
        PopFit(i) = objFun(PopPos(i,:));
        theta = (1 - 10*i/nPop) * pi;
        r = alpha * exp(beta * theta / 3);
        x(i) = r * cos(theta);
        y(i) = r * sin(theta);

        if PopFit(i) <= bestFitness
            bestFitness  = PopFit(i);
            bestPosition = PopPos(i,:);
        end
    end

    [~, windex] = max(PopFit);
    WorstX = PopPos(windex,:);
    Dis = abs(mean(mean((PopPos - bestPosition) ./ (WorstX - bestPosition + eps))));
    Lx  = abs(randn) * rand;

    convergenceCurve = zeros(maxItr,1);
    s = 0; z = 0; t = 0;

    %% Main Loop
    for It = 1:maxItr

        C = (1 - It/maxItr);
        T = (1 - (sin((pi*It)/(2*maxItr))))^(It/maxItr);

        if t < 15
            die = 0.02 * T;
        else
            die = 0.02;
            C = 0.8;
        end

        for i = 1:nPop

            Q = 0;
            F = sign(0.5 - rand);
            E = 1*T + rand;
            R1 = rand(1,dim);
            R4 = rand(1,dim);
            r1 = rand;
            r2 = rand;
            S = sin(pi * R4 * C);

            k = randperm(nPop);
            PopPosrand = PopPos(k,:);
            PopFitrand = PopFit(k);

            if rand > C

                J = abs(mean((PopPos(i,:) - bestPosition) ./ (WorstX - bestPosition + eps)));

                if rand > C
                    Cy = 1/(pi*(1+C^2));
                    if J > Dis
                        newPos = bestPosition + F*S.*(bestPosition - PopPos(i,:));
                    else
                        if Dis*Lx < J
                            newPos = bestPosition*(1+T^5*Cy*E) + F*(S.*(bestPosition - PopPos(i,:)));
                        else
                            newPos = bestPosition*(1+T^5*normrnd(0,C^2)) + F*(S.*(bestPosition - PopPos(i,:)));
                        end
                    end
                else
                    if rand > C
                        if mod(i,2)==1
                            r3 = rand;
                            step = (bestPosition - E*PopPos(i,:));
                            newPos = C/It*(r1*bestPosition - r3*PopPos(i,:)) ...
                                + T^2*levyFlight(1,dim).*abs(step);
                        else
                            R2 = rand(1,dim);
                            R3 = rand(1,dim);
                            step = PopPos(i,:) - E*bestPosition;
                            DE = C*F;
                            newPos = 0.5*(bestPosition + PopPosrand(1,:)) ...
                                + DE*(2*R1.*(step) - R2/2.*(DE*R3 - 1));
                        end
                    else
                        if rand < rand
                            if J < Dis
                                V = 2*(rand*(mean(PopPos)-PopPos(i,:)) ...
                                    + rand*(bestPosition-PopPos(i,:)));
                            else
                                V = 2*(rand*(PopPosrand(2,:)-PopPosrand(3,:)) ...
                                    + rand*(PopPosrand(1,:)-PopPos(i,:)));
                            end

                            if PopFit(i) <= PopFitrand(i)
                                step = PopPos(i,:) - E*PopPosrand(i,:);
                                base = PopPos(i,:);
                            else
                                step = PopPosrand(i,:) - E*PopPos(i,:);
                                base = PopPosrand(i,:);
                            end

                            if mod(i,2)==1
                                newPos = (base + T^2*y(i)*(1-R1).*abs(step)) ...
                                    + F*R1.*step/2 + V*J/It;
                            else
                                newPos = (base + T^2*x(i)*(1-R1).*abs(step)) ...
                                    + F*R1.*step/2 + V*J/It;
                            end

                            s = s + 1;
                            if s > 10
                                lesp1 = r1*PopPos(randperm(nPop,1),:) ...
                                    + (1-r1)*PopPos(randperm(nPop,1),:);
                                newPos = round(lesp1) ...
                                    + F*r1*R1/(It^4).*newPos;
                                s = 0;
                            end
                        else
                            [~, index] = sort(PopFit);
                            A2 = randi(4);
                            A1 = randi(4);
                            D = [PopPos(index(1:3),:); mean(PopPos)];
                            B = D(A1,:);
                            Rt1 = randperm(360,dim)*pi/360;
                            Rt2 = randperm(360,dim)*pi/360;
                            w = 1-((exp(It/maxItr)-1)/(exp(1)-1))^2;

                            if rand < 0.33
                                newPos = B + 2*w*F*(cos(Rt1)).*(sin(Rt2)).*(B-PopPos(i,:));
                            elseif rand < 0.66
                                newPos = B + 2*w*F*(sin(Rt1)).*(cos(Rt2)).*(B-PopPos(i,:));
                            else
                                newPos = B + 2*w*F*(cos(Rt2)).*(B-PopPos(i,:));
                            end

                            if A2 == 4, Q = 1; end

                            z = z + 1;
                            if z > 5
                                newPos = bestPosition.* ...
                                    (1-(1-1./(PopPos(randperm(nPop,1),:)+eps)).*R1);
                                z = 0;
                            end
                        end
                    end
                end

            else

                if rand > C
                    if rand > C
                        newPos = PopPosrand(3,:) + abs(randn)* ...
                            (bestPosition - PopPos(i,:) ...
                            + PopPosrand(1,:) - PopPosrand(2,:));
                    else
                        Z2 = rand(1,dim) < rand;
                        newPos = Z2.*(PopPosrand(3,:) ...
                            + abs(randn)*(PopPosrand(1,:) - PopPosrand(2,:))) ...
                            + (1-Z2).*PopPos(i,:);
                    end
                else
                    Z1 = rand < rand;
                    newPos = PopPos(i,:) ...
                        + (Z1*abs(randn)*((bestPosition+PopPosrand(1,:))/2 ...
                        - PopPosrand(2,:)) ...
                        + rand/2*(PopPosrand(3,:) - PopPosrand(4,:)));
                end

                if rand > C || t > 0.8*nPop
                    for j = 1:dim
                        if rand >= (0.2*C+0.2)
                            newPos(j) = PopPos(i,j);
                        end
                    end
                end
            end

            if rand < die
                if rand > C
                    newPos = rand(1,dim).*(ub-lb)+lb;
                else
                    best = bestPosition*(levyFlight(1,1)*(r1>r2) ...
                        + abs(randn)*(r1<=r2));
                    Upc = max(best);
                    Lowc = min(best);
                    newPos = rand(1,dim)*(Upc-Lowc)+Lowc;
                end
            end

            newPos = min(max(newPos,lb),ub);
            newFit = objFun(newPos);

            if newFit < PopFit(i)
                PopFit(i) = newFit;
                PopPos(i,:) = newPos;

                if Q == 1
                    [~, index] = sort(PopFit);
                    PopPos(index(end),:) = PopPos(i,:);
                    PopFit(index(end)) = PopFit(i);
                end
                t = 0;
            else
                t = t + 1;
            end

            if PopFit(i) <= bestFitness
                bestFitness  = PopFit(i);
                bestPosition = PopPos(i,:);
            end

            [~, windex] = max(PopFit);
            WorstX = PopPos(windex,:);
        end

        convergenceCurve(It) = bestFitness;
    end

end

%% Lévy Flight
function L = levyFlight(n,m)
    Beta = 1.5;
    sigma = (gamma(1+Beta)*sin(pi*Beta/2) / ...
        (gamma((1+Beta)/2)*Beta*2^((Beta-1)/2)))^(1/Beta);
    u = randn(n,m) * sigma;
    v = randn(n,m);
    L = 0.05 * u ./ abs(v).^(-Beta);
end
