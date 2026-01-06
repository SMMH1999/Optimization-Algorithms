function [bestScore, bestPos, curve] = SSA(LB, UB, Dim, Salp_no, Max_iter, objective)
    %% Initialize Parameters
    if size(UB, 2) == 1
        UB = ones(1, Dim) * UB;
        LB = ones(1, Dim) * LB;
    end
    curve = zeros(Max_iter, 1);

    %% Initialize positions and fitness of salps
    SalpPositions = Population_Generator(Salp_no, Dim, UB, LB);
    SalpFitness = inf(1, Salp_no);
    bestPos = zeros(1, Dim);
    bestScore = inf;

    %% Calculate the objective function for each salp
    for i = 1:size(SalpPositions, 1)
        SalpFitness(i) = objective(SalpPositions(i, :));
    end

    %% Sorting the Population
    [sorted_salps_fitness, sorted_indexes] = sort(SalpFitness);
    for newindex = 1:Salp_no
        Sorted_salps(newindex, :) = SalpPositions(sorted_indexes(newindex), :);
    end
    bestPos = Sorted_salps(1, :);
    bestScore = sorted_salps_fitness(1);

    %% Main loop
    % start from the second iteration since the first iteration was dedicated to calculating the fitness of salps
    l = 2;
    while l < Max_iter + 1
        c1 = 2 * exp(-(4 * l / Max_iter) ^ 2);
        for i = 1:size(SalpPositions, 1)
            SalpPositions= SalpPositions';
            if i <= Salp_no / 2
                for j = 1:1:Dim
                    if rand() < 0.5
                        SalpPositions(j, i) = bestPos(j) + c1 * ((UB(j) - LB(j)) * rand() + LB(j));
                    else
                        SalpPositions(j, i) = bestPos(j) - c1 * ((UB(j) - LB(j)) * rand() + LB(j));
                    end
                end
            elseif i > Salp_no / 2 && i < Salp_no + 1
                point1 = SalpPositions(:, i-1);
                point2 = SalpPositions(:, i);
                SalpPositions(:, i) = (point2 + point1) / 2;
            end
            SalpPositions = SalpPositions';
        end
        %% Return back the search agents that go beyond the boundaries of the search space
        SalpPositions = min(max(SalpPositions, LB), UB);

        %% Calculate the objective function for each salp
        for i = 1:size(SalpPositions, 1)
            SalpFitness(i) = objective(SalpPositions(i, :));
        end

        for i = 1:size(SalpPositions, 1)
            if SalpFitness(i) < bestScore
                bestPos = SalpPositions(i, :);
                bestScore = SalpFitness(i);
            end
        end
        curve(l) = bestScore;
        l = l + 1;
    end
end
