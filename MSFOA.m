function [xposbest, fvalbest, Curve] = MSFOA(Npop, Max_it, lb, ub, nD, fobj)
% Modified Starfish Optimization Algorithm (MSFOA)
% Enhancements added to the original SFOA:
% 1. Adaptive step size
% 2. Momentum vector
% 3. Multi-candidate local search
% 4. Dynamic exploration–exploitation balance
%
% NOTE:
% This implementation is inspired by and should be cited as:
% H. Akbulut, "A modified starfish optimization algorithm (M-SFOA) for global
% optimization problems and its application to heart disease risk prediction,"
% Expert Systems with Applications, vol. 307, 2026, Article 131088.
% https://doi.org/10.1016/j.eswa.2026.131088

%% Initial Parameters
GP = 0.5;     % Base parameter

% Boundary handling
if size(ub, 2) == 1
    lb = lb * ones(1, nD); 
    ub = ub * ones(1, nD);
end

% Enhanced parameters
step_size_initial = 1.0;      % Initial step size
step_size_final   = 0.01;     % Final step size
momentum_factor   = 0.7;      % Momentum factor
prev_best = zeros(1, nD);     % Previous best solution

%% Initial Population
fvalbest = inf;
Curve = zeros(1, Max_it);
Xpos = rand(Npop, nD) .* (ub - lb) + lb;

% Fitness evaluation
Fitness = zeros(1, Npop);
for i = 1:Npop
    Fitness(i) = feval(fobj, Xpos(i, :));
end

[fvalbest, order] = min(Fitness);    % Global best fitness
xposbest = Xpos(order, :);           % Global best position
prev_best = xposbest;                % Store previous best solution

newX = zeros(Npop, nD);

%% Optimization Loop
for T = 1:Max_it
    % Adaptive step size computation (exponential decay)
    current_step_size = step_size_initial * ...
        (step_size_final / step_size_initial)^(T / Max_it);
    
    % Momentum vector computation
    momentum = xposbest - prev_best;
    prev_best = xposbest;
    
    theta = pi/2 * T / Max_it;
    tEO = (Max_it - T)/Max_it * cos(theta);
    
    if rand < GP  % Exploration phase (starfish movement)
        for i = 1:Npop
            if nD > 5
                % For dimensionality greater than 5 (momentum-enhanced version)
                jp1 = randperm(nD, 5);
                for j = 1:5
                    pm = (2*rand - 1) * pi;
                    if rand < GP
                        % Momentum-enhanced exploratory movement
                        newX(i, jp1(j)) = Xpos(i, jp1(j)) + current_step_size * ...
                            (pm * (xposbest(jp1(j)) - Xpos(i, jp1(j))) * cos(theta) + ...
                             momentum_factor * momentum(jp1(j)));
                    else
                        newX(i, jp1(j)) = Xpos(i, jp1(j)) - pm * ...
                            (xposbest(jp1(j)) - Xpos(i, jp1(j))) * sin(theta);
                    end
                    % Boundary check
                    if newX(i, jp1(j)) > ub(jp1(j)) || newX(i, jp1(j)) < lb(jp1(j))
                        newX(i, jp1(j)) = Xpos(i, jp1(j));
                    end
                end
            else
                % For dimensionality 5 or less (multi-candidate approach)
                jp2 = ceil(nD * rand);
                im = randperm(Npop);
                rand1 = 2*rand - 1;
                rand2 = 2*rand - 1;
                
                % Primary candidate
                candidate1 = tEO * Xpos(i, jp2) + ...
                    rand1 * (Xpos(im(1), jp2) - Xpos(i, jp2)) + ...
                    rand2 * (Xpos(im(2), jp2) - Xpos(i, jp2));
                
                % Secondary candidate (with momentum)
                candidate2 = Xpos(i, jp2) + current_step_size * ...
                    (rand1 * (xposbest(jp2) - Xpos(i, jp2)) + ...
                     momentum_factor * momentum(jp2));
                
                % Third candidate (random)
                candidate3 = Xpos(i, jp2) + current_step_size * randn;
                
                % Select the best candidate
                candidates = [candidate1, candidate2, candidate3];
                [~, best_idx] = min(abs(candidates - xposbest(jp2)));
                newX(i, jp2) = candidates(best_idx);
                
                % Boundary check
                if newX(i, jp2) > ub(jp2) || newX(i, jp2) < lb(jp2)
                    newX(i, jp2) = Xpos(i, jp2);
                end
            end
            newX(i, :) = max(min(newX(i, :), ub), lb);  % Boundary enforcement
        end
    else  % Exploitation phase (starfish feeding behavior)
        df = randperm(Npop, 5);
        dm = zeros(5, nD);
        for k = 1:5
            dm(k, :) = xposbest - Xpos(df(k), :);
        end
        
        for i = 1:Npop
            r1 = rand; 
            r2 = rand;
            kp = randperm(5, 2);
            
            % Primary movement vector
            move1 = Xpos(i, :) + r1 .* dm(kp(1), :) + r2 .* dm(kp(2), :);
            
            % Momentum-enhanced secondary vector
            move2 = Xpos(i, :) + current_step_size * ...
                (r1 .* dm(kp(1), :) + momentum_factor * momentum);
            
            % Multi-candidate evaluation
            if i == Npop
                % Regeneration movement (with adaptive step size)
                move3 = exp(-T * Npop / Max_it) .* Xpos(i, :);
            else
                % Random perturbation
                move3 = Xpos(i, :) + current_step_size * randn(1, nD);
            end
            
            % Evaluate candidates and select the best one
            candidates = [move1; move2; move3];
            candidate_fitness = zeros(1, 3);
            for c = 1:3
                candidate_fitness(c) = feval(fobj, candidates(c, :));
            end
            [~, best_candidate] = min(candidate_fitness);
            newX(i, :) = candidates(best_candidate, :);
            
            % Boundary check
            newX(i, :) = max(min(newX(i, :), ub), lb);
        end
    end
    
    % Fitness evaluation and update
    for i = 1:Npop
        newFit = feval(fobj, newX(i, :));
        if newFit < Fitness(i)
            Fitness(i) = newFit;
            Xpos(i, :) = newX(i, :);
            if newFit < fvalbest
                fvalbest = newFit;
                xposbest = Xpos(i, :);
            end
        end
    end
    
    % Update convergence curve
    Curve(T) = fvalbest;
    
    % Display progress
    if mod(T, 50) == 0 || T == 1
        fprintf('Iteration %d, Best Value: %f\n', T, fvalbest);
    end
end
end
