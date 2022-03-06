classdef OperatorDeroa < OperatorJADE
    methods (Static)
        function [offspring, r1_task, calls] = generate(callfun, population,p, Tasks, k, rmp)
            if isempty(population{k})
                offspring = population{k};
                calls = 0;
                return;
            end
            Individual_class = class(population{k}(1));
            
            % get top 100p% individuals
            factorial_costs = [];
            for i = 1:length(population{k})
                factorial_costs(i) = population{k}(i).factorial_costs;
            end
            [~, rank] = sort(factorial_costs);
            pop_pbest = rank(1:round(p * length(population{k})));
            
            r1_task = zeros(1, length(population{k}));
            for i = 1:length(population{k})
                offspring(i) = feval(Individual_class);

                pbest = pop_pbest(randi(length(pop_pbest)));
                
                rnd = randperm(length(population{k}), 3);
                x1 = rnd(1); x2 = rnd(2);x3 = rnd(3);

                r = rand();
                for t = 1:length(Tasks)
                    if r <= sum(rmp(k, 1:t))
                        r1_task(i) = t;
                        break;
                    end
                end

                offspring(i) = OperatorJADE.mutate_current_pbest_1(offspring(i), population{r1_task(i)}(x1), population{r1_task(i)}(pbest),population{k}(x2), population{k}(x3));
                offspring(i) = OperatorJADE.crossover(offspring(i), population{k}(i));

                vio_low = find(offspring(i).rnvec < 0);
                offspring(i).rnvec(vio_low) = (population{k}(i).rnvec(vio_low) + 0) / 2;
                vio_up = find(offspring(i).rnvec > 1);
                offspring(i).rnvec(vio_up) = (population{k}(i).rnvec(vio_up) + 1) / 2;
            end
            if callfun
                [offspring, calls] = evaluate(offspring, Tasks(k), 1);
            else
                calls = 0;
            end
        end
    end
end
