function [bestFitness, bestPosition, convergenceCurve] = AVOA(lb, ub, dim, nPop, maxItr, objFun)

    %% -------------------- Initialization -------------------- %%
    % Ensure row vectors
    lb = lb(:)';
    ub = ub(:)';

    % Scalar bound support
    if isscalar(lb)
        lb = repmat(lb, 1, dim);
        ub = repmat(ub, 1, dim);
    end

    % Initialize population
    X = rand(nPop, dim) .* (ub - lb) + lb;

    % Best solutions
    bestPosition  = zeros(1, dim);
    bestFitness   = inf;
    secondBestPos = zeros(1, dim);
    secondBestFit = inf;

    convergenceCurve = zeros(1, maxItr);

    %% -------------------- Control Parameters -------------------- %%
    p1 = 0.6;
    p2 = 0.4;
    p3 = 0.6;
    alpha = 0.8;
    beta  = 0.2;
    gammaVal = 2.5;

    %% ==================== Main Loop ==================== %%
    for t = 1:maxItr

        %% -------- Fitness Evaluation -------- %%
        for i = 1:nPop
            fitness = objFun(X(i,:));

            if fitness < bestFitness
                secondBestFit = bestFitness;
                secondBestPos = bestPosition;

                bestFitness = fitness;
                bestPosition = X(i,:);

            elseif fitness > bestFitness && fitness < secondBestFit
                secondBestFit = fitness;
                secondBestPos = X(i,:);
            end
        end

        %% -------- Adaptive Control Parameter -------- %%
        a = (rand()*4 - 2) * ...
            ((sin((pi/2)*(t/maxItr))^gammaVal) + ...
            cos((pi/2)*(t/maxItr)) - 1);

        P1 = (2*rand + 1)*(1 - t/maxItr) + a;

        %% -------- Position Update -------- %%
        for i = 1:nPop

            F = P1 * (2*rand - 1);

            % Random vulture selection (roulette)
            if rand <= alpha
                randomVulture = bestPosition;
            else
                randomVulture = secondBestPos;
            end

            Xi = X(i,:);

            % ================= Exploration =================
            if abs(F) >= 1

                if rand < p1
                    Xi = randomVulture - ...
                        abs((2*rand)*randomVulture - Xi) * F;
                else
                    Xi = randomVulture - F + ...
                        rand*( (ub-lb).*rand(1,dim) + lb );
                end

                % ================= Exploitation =================
            else

                % ---- Phase 1 ----
                if abs(F) < 0.5

                    if rand < p2
                        A = bestPosition - ...
                            ((bestPosition.*Xi) ./ ...
                            (bestPosition - Xi.^2)) * F;

                        B = secondBestPos - ...
                            ((secondBestPos.*Xi) ./ ...
                            (secondBestPos - Xi.^2)) * F;

                        Xi = (A + B) / 2;

                    else
                        levy = levyFlight(dim);
                        Xi = randomVulture - ...
                            abs(randomVulture - Xi) * F .* levy;
                    end

                    % ---- Phase 2 ----
                else

                    if rand < p3
                        Xi = abs((2*rand)*randomVulture - Xi) ...
                            *(F + rand) ...
                            - (randomVulture - Xi);
                    else
                        s1 = randomVulture .* ...
                            (rand()*Xi/(2*pi)) .* cos(Xi);
                        s2 = randomVulture .* ...
                            (rand()*Xi/(2*pi)) .* sin(Xi);
                        Xi = randomVulture - (s1 + s2);
                    end

                end
            end

            X(i,:) = Xi;
        end

        %% -------- Boundary Control -------- %%
        X = min(max(X, lb), ub);

        %% -------- Convergence Record -------- %%
        convergenceCurve(t) = bestFitness;

    end

end

%% ================= Local Function ================= %%
function step = levyFlight(d)

    beta = 3/2;

    sigma = (gamma(1+beta)*sin(pi*beta/2) / ...
        (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

    u = randn(1,d) * sigma;
    v = randn(1,d);

    step = u ./ abs(v).^(1/beta);

end
