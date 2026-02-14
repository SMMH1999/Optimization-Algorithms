function [bestFitness, bestPosition, convergenceCurve] = TSAO(lb, ub, dim, nPop, maxItr, objFun)
    % Transit Search (TS)
    % -------------------------------------------------------------------------
    % Transit Search (TS) is a physics-inspired metaheuristic optimization
    % algorithm based on galactic structures, stellar luminosity, planetary
    % transit detection, and signal amplification mechanisms.
    %
    % INPUTS:
    %   lb        : [1×dim] or scalar, lower bound(s) of decision variables
    %   ub        : [1×dim] or scalar, upper bound(s) of decision variables
    %   dim       : Scalar, number of decision variables
    %   nPop      : Scalar, number of stars (population size) (ns)
    %   maxItr    : Scalar, maximum number of iterations (maxcycle)
    %   objFun    : Function handle, objective function (minimization)
    %
    % OUTPUTS:
    %   bestFitness      : Scalar, best objective value found
    %   bestPosition     : [1×dim] vector, best solution found
    %   convergenceCurve : [1×maxItr] vector, best fitness per iteration
    %
    % Tunable Parameters (internal):
    %   SN  : Number of neighboring samples per star (fixed = 5)
    %
    % -------------------------------------------------------------------------

    %% Parameters
    ns = nPop;
    SN = 5;

    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    end

    %% Initialization

    Empty.Location = [];
    Empty.Cost = inf;

    Galaxy_Center = Empty;
    region = repmat(Empty, ns*SN, 1);
    selected_regions = repmat(Empty, ns, 1);
    Stars = repmat(Empty, ns, 1);
    Stars2 = repmat(Empty, ns, 1);
    Best_Planets = repmat(Empty, ns, 1);
    SN_P = repmat(Empty, SN, 1);

    Stars_sorted = zeros(ns,1);
    Ranks = 1:ns;
    Stars_Ranks = zeros(ns,1);
    Star_RanksNormal = zeros(ns,1);
    Distance = zeros(ns,1);
    Luminosity = zeros(ns,1);
    Luminosity_new = zeros(ns,1);
    Transit0 = zeros(ns,1);

    convergenceCurve = zeros(1, maxItr);

    %% Galaxy Phase

    Galaxy_Center.Location = lb + rand(1,dim).*(ub-lb);
    Galaxy_Center.Cost = objFun(Galaxy_Center.Location);

    for l = 1:(ns*SN)
        if randi(2)==1
            difference = rand.*Galaxy_Center.Location - (lb + rand(1,dim).*(ub-lb));
        else
            difference = rand.*Galaxy_Center.Location + (lb + rand(1,dim).*(ub-lb));
        end
        Noise = (rand(1,dim).^3).*(lb + rand(1,dim).*(ub-lb));
        region(l).Location = Galaxy_Center.Location + difference - Noise;
        region(l).Location = max(region(l).Location, lb);
        region(l).Location = min(region(l).Location, ub);
        region(l).Cost = objFun(region(l).Location);
    end

    [~, index] = sort([region.Cost]);
    for i = 1:ns
        selected_regions(i) = region(index(i));
        for k = 1:SN
            if randi(2)==1
                difference = rand.*selected_regions(i).Location - rand.*(lb + rand(1,dim).*(ub-lb));
            else
                difference = rand.*selected_regions(i).Location + rand.*(lb + rand(1,dim).*(ub-lb));
            end
            Noise = (rand(1,dim).^3).*(lb + rand(1,dim).*(ub-lb));
            new.Location = selected_regions(i).Location + difference - Noise;
            new.Location = max(new.Location, lb);
            new.Location = min(new.Location, ub);
            new.Cost = objFun(new.Location);
            if new.Cost < Stars(i).Cost
                Stars(i) = new;
            end
        end
    end

    Best_Planets = Stars;
    [~, index] = sort([Best_Planets.Cost]);
    Best_Planet = Best_Planets(index(1));

    Telescope.Location = lb + rand(1,dim).*(ub-lb);

    for i = 1:ns
        Stars_sorted(i) = Stars(i).Cost;
    end
    Stars_sorted = sort(Stars_sorted);

    for i = 1:ns
        for ii = 1:ns
            if Stars(i).Cost == Stars_sorted(ii)
                Stars_Ranks(i) = Ranks(ii);
                Star_RanksNormal(i) = Stars_Ranks(i)/ns;
            end
        end
        Distance(i) = norm(Stars(i).Location - Telescope.Location);
        Luminosity(i) = Star_RanksNormal(i)/(Distance(i)^2 + eps);
    end
    Luminosity_new = Luminosity;

    %% Main Loop

    for it = 1:maxItr

        Transit = Transit0;
        Luminosity = Luminosity_new;

        %% Transit Phase
        for i = 1:ns
            difference = (2*rand-1).*Stars(i).Location;
            Noise = (rand(1,dim).^3).*(lb + rand(1,dim).*(ub-lb));
            Stars2(i).Location = Stars(i).Location + difference - Noise;
            Stars2(i).Location = max(Stars2(i).Location, lb);
            Stars2(i).Location = min(Stars2(i).Location, ub);
            Stars2(i).Cost = objFun(Stars2(i).Location);
        end

        for i = 1:ns
            Stars_sorted(i) = Stars2(i).Cost;
        end
        Stars_sorted = sort(Stars_sorted);

        for i = 1:ns
            for ii = 1:ns
                if Stars2(i).Cost == Stars_sorted(ii)
                    Stars_Ranks(i) = Ranks(ii);
                    Star_RanksNormal(i) = Stars_Ranks(i)/ns;
                end
            end
            Distance(i) = norm(Stars2(i).Location - Telescope.Location);
            Luminosity_new(i) = Star_RanksNormal(i)/(Distance(i)^2 + eps);
            if Luminosity_new(i) < Luminosity(i)
                Transit(i) = 1;
            end
        end
        Stars = Stars2;

        %% Location Phase (Exploration)
        for i = 1:ns
            if Transit(i)==1

                Luminosity_Ratio = Luminosity_new(i)/Luminosity(i);
                Planet.Location = (rand.*Telescope.Location + Luminosity_Ratio.*Stars(i).Location)/2;

                for k = 1:SN
                    zone = randi(3);
                    if zone==1
                        new.Location = Planet.Location - (2*rand-1).*(lb + rand(1,dim).*(ub-lb));
                    elseif zone==2
                        new.Location = Planet.Location + (2*rand-1).*(lb + rand(1,dim).*(ub-lb));
                    else
                        new.Location = Planet.Location + (2*rand(1,dim)-1).*(lb + rand(1,dim).*(ub-lb));
                    end
                    new.Location = max(new.Location, lb);
                    new.Location = min(new.Location, ub);
                    SN_P(k) = new;
                end

                SUM = zeros(1,dim);
                for k = 1:SN
                    SUM = SUM + SN_P(k).Location;
                end
                new.Location = SUM/SN;
                new.Cost = objFun(new.Location);

                if new.Cost < Best_Planets(i).Cost
                    Best_Planets(i) = new;
                end

            else

                Neighbor.Location = (rand.*Stars(i).Location + rand.*(lb + rand(1,dim).*(ub-lb)))/2;

                for k = 1:SN
                    zone = randi(3);
                    if zone==1
                        Neighbor.Location = Neighbor.Location - (2*rand-1).*(lb + rand(1,dim).*(ub-lb));
                    elseif zone==2
                        Neighbor.Location = Neighbor.Location + (2*rand-1).*(lb + rand(1,dim).*(ub-lb));
                    else
                        Neighbor.Location = Neighbor.Location + (2*rand(1,dim)-1).*(lb + rand(1,dim).*(ub-lb));
                    end
                    Neighbor.Location = max(Neighbor.Location, lb);
                    Neighbor.Location = min(Neighbor.Location, ub);
                    Neighbor.Cost = objFun(Neighbor.Location);
                    SN_P(k) = Neighbor;
                end

                SUM = zeros(1,dim);
                for k = 1:SN
                    SUM = SUM + SN_P(k).Location;
                end
                Neighbor.Location = SUM/SN;
                Neighbor.Cost = objFun(Neighbor.Location);

                if Neighbor.Cost < Best_Planets(i).Cost
                    Best_Planets(i) = Neighbor;
                end
            end
        end

        %% Signal Amplification (Exploitation)
        for i = 1:ns
            for k = 1:SN
                if randi(2)==1
                    Power = randi(SN*ns);
                    Coefficient = 2*rand;
                    Noise = (rand(1,dim).^Power).*(lb + rand(1,dim).*(ub-lb));
                else
                    Power = randi(SN*ns);
                    Coefficient = 2*rand;
                    Noise = -(rand(1,dim).^Power).*(lb + rand(1,dim).*(ub-lb));
                end

                if randi(2)==1
                    new.Location = Best_Planets(i).Location - Coefficient.*Noise;
                else
                    new.Location = (rand.*Best_Planets(i).Location) - Coefficient.*Noise;
                end

                new.Location = max(new.Location, lb);
                new.Location = min(new.Location, ub);
                new.Cost = objFun(new.Location);

                if new.Cost < Best_Planets(i).Cost
                    Best_Planets(i) = new;
                end
            end

            if Best_Planets(i).Cost < Best_Planet.Cost
                Best_Planet = Best_Planets(i);
            end
        end

        bestFitness = Best_Planet.Cost;
        bestPosition = Best_Planet.Location;
        convergenceCurve(it) = bestFitness;

    end

end
