classdef Deroa2 < Algorithm
    % <Multi> <None>

    properties (SetAccess = private)
        p=0.1
        H=50
        alpha = 0.5
        beta = 0.2
        gama = 0.85
        rmp0 = 0.3
    end

    methods
        function parameter = getParameter(obj)
            parameter = {'p: 100p% top as pbest', num2str(obj.p), ...
                        'H: success memory size', num2str(obj.H), ...
                        'alpha', num2str(obj.alpha), ...
                        'beta', num2str(obj.beta), ...
                        'gama', num2str(obj.gama), ...
                        'rmp0: Initial random mating probability', num2str(obj.rmp0)};
        end

        function obj = setParameter(obj, parameter_cell)
            count = 1;
            obj.p = str2double(parameter_cell{count}); count = count + 1;
            obj.H = str2double(parameter_cell{count}); count = count + 1;
            obj.alpha = str2double(parameter_cell{count}); count = count + 1;
            obj.beta = str2double(parameter_cell{count}); count = count + 1;
            obj.gama = str2double(parameter_cell{count}); count = count + 1;
            obj.rmp0 = str2double(parameter_cell{count}); count = count + 1;
        end

        function data = run(obj, Tasks, run_parameter_list)
            sub_pop = run_parameter_list(1);
            iter_num = run_parameter_list(2);
            eva_num = run_parameter_list(3) * length(Tasks);
            tic

            pop_size = sub_pop * length(Tasks);
            population = {};
            fnceval_calls = 0;
%             HR = []; % HR is used to store the historical rewards
            if iter_num == inf
                gen = (eva_num - (pop_size * length(Tasks) - 1)) / pop_size;
                T = obj.beta * gen;
                delta_rmp = 1 / gen;
            else
                T = obj.beta * iter_num;
                delta_rmp = 1 / iter_num;
            end
            rmp = obj.rmp0 * ones(length(Tasks), length(Tasks)) / (length(Tasks) - 1);
            rmp(logical(eye(size(rmp)))) = (1 - obj.rmp0);

            for t = 1:length(Tasks)
                for i = 1:sub_pop
                    population{t}(i) = IndividualJADE();
                    population{t}(i).rnvec = rand(1, max([Tasks.dims]));
                end
                [population{t}, calls] = evaluate(population{t}, Tasks(t), 1);
                fnceval_calls = fnceval_calls + calls;

                % initialize parameter
                H_idx(t) = 1;
                MF{t} = 0.5 .* ones(1, obj.H);
                MCR{t} = 0.5 * ones(1, obj.H);
               
                
                [bestobj(t), idx] = min([population{t}.factorial_costs]);
                data.bestX{t} = population{t}(idx).rnvec;
                data.convergence(t, 1) = bestobj(t);
            end

            generation = 1;
            while generation < iter_num && fnceval_calls < eva_num
                generation = generation + 1;
                
                for t = 1:length(Tasks)               
                   % calculate individual F and CR
                    for i = 1:length(population{t})
                        idx = randi(obj.H);
                        uF = MF{t}(idx);
                        population{t}(i).F = uF + 0.1 * tan(pi * (rand - 0.5));
                        while (population{t}(i).F <= 0)
                            population{t}(i).F = uF + 0.1 * tan(pi * (rand - 0.5));
                        end
                        population{t}(i).F(population{t}(i).F > 1) = 1;

                        uCR = MCR{t}(idx);
                        population{t}(i).CR = normrnd(uCR, 0.1);
                        population{t}(i).CR(population{t}(i).CR > 1) = 1;
                        population{t}(i).CR(population{t}(i).CR < 0) = 0;
                    end
                    
%                 % Select the k-th task to optimize
%                 if generation <= T
%                     k = unidrnd(length(Tasks));
%                 else
%                     weights = obj.gama.^(generation - 3:-1:0);
%                     sum_weights = sum(weights);
%                     for t = 1:length(Tasks)
%                         mean_R(t) = sum(weights .* HR(t, :)) / sum_weights;
%                     end
%                     % The selection probability
%                     prob(generation, :) = obj.prob_min / length(Tasks) + (1 - obj.prob_min) * mean_R ./ (sum(mean_R));
%                     % Determine the a task based on the selection probability using roulette wheel method
%                     r = rand;
%                     for t = 1:length(Tasks)
%                         if r <= sum(prob(generation, 1:t))
%                             k = t;
%                             break;
%                         end
%                     end
%                 end

                     % generate
                     [offspring, r1_task, calls] = OperatorDeroa2.generate(1,population, obj.p,Tasks,rmp,t);
                     fnceval_calls = fnceval_calls + calls;

                     % selection
                     fit_old = [population{t}.factorial_costs];
                     replace = [population{t}.factorial_costs] > [offspring.factorial_costs];
                     population{t}(replace) = offspring(replace);
                     fit_new = [population{t}.factorial_costs];
                
                     % update archive
                    arc{t} = [arc{t}, population{t}(replace)];
                    if length(arc{t}) > length(population{t})
                        rnd = randperm(length(arc{t}));
                        arc{t} = arc{t}(rnd(1:length(population{t})));
                    end
              
                    % calculate SF SCR
                    SF = [population{t}(replace).F];
                    SCR = [population{t}(replace).CR];
                    dif = abs([population{t}(replace).factorial_costs] - [offspring(replace).factorial_costs]);
                    dif = dif ./ sum(dif);

                    % update MF MCR
                    if ~isempty(SF)
                      MF{t}(H_idx(t)) = (dif * (SF'.^2)) / (dif * SF');
                      MCR{t}(H_idx(t)) = (dif * (SCR'.^2)) / (dif * SCR');
                    else
                      MF{t}(H_idx(t)) = MF{t}(mod(H_idx(t) + obj.H - 2, obj.H) + 1);
                      MCR{t}(H_idx(t)) = MCR{t}(mod(H_idx(t) + obj.H - 2, obj.H) + 1);
                    end
                    H_idx(t) = mod(H_idx(t), obj.H) + 1;
              
                    % calculate the reward
                   R_p = max((fit_old - fit_new) ./ (fit_old), 0);
                   R_b = max((min(bestobj) - min(fit_new)) / (min(bestobj)), 0);
                   R = zeros(length(Tasks), 1);
                   for k = 1:length(Tasks)
                       if k == t %The main task
                          R(k) = obj.alpha * R_b + (1 - obj.alpha) * (sum(R_p) / pop_size);
                       else % The auxiliary task
                          index = find(r1_task == k);
                          if isempty(index)
                            R(k) = 0;
                          else
                            [~, minid] = min(fit_new);
                            R(k) = obj.alpha * (r1_task(minid) == k) * R_b + (1 - obj.alpha) * (sum(R_p(index)) / length(index));
                          end
                       end
                   end
%                    HR = [HR, R];
                   
                   % update rmp
                   for k = 1:length(Tasks)
                       if k == t
                          continue;
                       end
                       if R(k) >= R(t)
                           rmp(t, k) = min(rmp(t, k) + delta_rmp, 1);
                           rmp(t, t) = max(rmp(t, t) - delta_rmp, 0);
                       else
                           rmp(t, k) = max(rmp(t, k) - delta_rmp, 0);
                           rmp(t, t) = min(rmp(t, t) + delta_rmp, 1);
                       end
                   end

                population{t}(replace) = offspring(replace);
                [bestobj_now, idx] = min([population{t}.factorial_costs]);
                if bestobj_now < bestobj(t)
                    bestobj(t) = bestobj_now;
                    data.bestX{t} = population{k}(idx).rnvec;
                end
                data.convergence(t, generation) = bestobj(t);
               end
            end
            data.bestX = uni2real(data.bestX, Tasks);
            data.clock_time = toc;
        end
    end
end
