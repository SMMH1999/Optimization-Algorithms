function [bestScore, bestPos, curve] = PSO(LB, UB, Dim, Particles_no, Max_iter, objective)
    %% Parameters
    w = 1;  % Inertia weight
    c1 = 2;  % Cognitive (personal) acceleration coefficient
    c2 = 2;  % Social (global) acceleration coefficient
    curve = zeros(Max_iter, 1);
    fitness = inf(1, Particles_no);

    %% Initialize the particles
    Positions = Population_Generator(Particles_no, Dim, UB, LB);
    velocities = zeros(Particles_no, Dim);  % Initial velocities
    personalBest_Positions = Positions;  % Best Positions of each particle
    personalBest_Fitnesss = inf(1, Particles_no);  % Best fitness of each particle
    bestScore = inf;  % Best score of the swarm
    bestPos = zeros(1, Dim);  % Best position of the swarm

    %% Main loop
    for Itr = 1 : Max_iter
        %% Calculate the objective function for each particle
        for i = 1 : Particles_no
            fitness(i) = objective(Positions(i,:));
        end

        %% Update personal best Positions and fitness
        betterIndices = fitness < personalBest_Fitnesss;
        personalBest_Fitnesss(betterIndices) = fitness(betterIndices);
        personalBest_Positions(betterIndices, :) = Positions(betterIndices, :);

        %% Update global best position and score
        [minFitness, minIndex] = min(fitness);
        if minFitness < bestScore
            bestScore = minFitness;
            bestPos = Positions(minIndex, :);
        end

        %% Update curve
        curve(Itr) = bestScore;
        if Itr > 1
            if bestScore > curve(Itr - 1)
                curve(Itr) = curve(Itr - 1);
            end
        end

        %% Update velocities and Positions
        for i = 1:Particles_no
            r1 = rand(1, Dim);
            r2 = rand(1, Dim);
            cognitiveComponent = c1 * r1 .* (personalBest_Positions(i, :) - Positions(i, :));
            socialComponent = c2 * r2 .* (bestPos - Positions(i, :));
            velocities(i, :) = w * velocities(i, :) + cognitiveComponent + socialComponent;
            Positions(i, :) = Positions(i, :) + velocities(i, :);

            % Apply the bounds to the Positions
            Positions(i, :) = max(min(Positions(i, :), UB), LB);
        end
    end
end
