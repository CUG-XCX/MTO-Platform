classdef OperatorDeroa2 < OperatorJADE
    methods (Static)
        function [offspring, r1_task, calls] = generate(callfun,population,p, Tasks,rmp,t)
            if isempty(population{t})
                offspring = population{t};
                calls = 0;
                return;
            end
            Individual_class = class(population{t}(1));
            
            % get top 100p% individuals
            factorial_costs = [];
            for i = 1:length(population{t})
                factorial_costs(i) = population{t}(i).factorial_costs;
            end
            [~, rank] = sort(factorial_costs);
            pop_pbest = rank(1:round(p * length(population{t})));
            
            r1_task = zeros(1, length(population));
            for i = 1:length(population{t})
                offspring(i) = feval(Individual_class);

                pbest = pop_pbest(randi(length(pop_pbest)));
                
                rnd = randperm(length(population{t}), 3);
                x1 = rnd(1); x2 = rnd(2);x3 = rnd(3);

                r = rand();
                for k = 1:length(Tasks) 
                    if r <= sum(rmp(t, 1:k))
                        r1_task(i) = k;               
                        break;
                    end                  
                end
                if r1_task(i)==0
                    r1_task(i) = t;
                end 
                population{r1_task(i)}(x1).F   
                offspring(i) = OperatorJADE.mutate_current_pbest_1(offspring(i), population{r1_task(i)}(x1), population{r1_task(i)}(pbest),population{t}(x2), population{t}(x3));
                offspring(i) = OperatorJADE.crossover(offspring(i), population{t}(i));
                
                vio_low = find(offspring(i).rnvec < 0);
                offspring(i).rnvec(vio_low) = (population{t}(i).rnvec(vio_low) + 0) / 2;
                vio_up = find(offspring(i).rnvec > 1);
                offspring(i).rnvec(vio_up) = (population{t}(i).rnvec(vio_up) + 1) / 2;
            end
            if callfun
                [offspring, calls] = evaluate(offspring, Tasks(t), 1);
            else
                calls = 0;
            end
        end
    end
end
