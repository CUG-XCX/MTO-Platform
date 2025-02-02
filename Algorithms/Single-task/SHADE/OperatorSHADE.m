classdef OperatorSHADE < OperatorJADE
    methods (Static)
        function [offspring, calls] = generate(callfun, Task, population, union, p)
            if length(population) <= 3
                offspring = population;
                calls = 0;
                return;
            end
            Individual_class = class(population(1));

            % get top 100p% individuals
            for i = 1:length(population)
                factorial_costs(i) = population(i).factorial_costs;
            end
            [~, rank] = sort(factorial_costs);
            pop_pbest = rank(1:round(p * length(population)));

            for i = 1:length(population)
                offspring(i) = feval(Individual_class);

                pbest = pop_pbest(randi(length(pop_pbest)));
                x1 = randi(length(population));
                while x1 == i
                    x1 = randi(length(population));
                end
                x2 = randi(length(union));

                offspring(i) = OperatorJADE.mutate_current_pbest_1(offspring(i), population(i), population(pbest), population(x1), union(x2));
                offspring(i) = OperatorJADE.crossover(offspring(i), population(i));

                offspring(i).rnvec(offspring(i).rnvec > 1) = 1;
                offspring(i).rnvec(offspring(i).rnvec < 0) = 0;
            end
            if callfun
                [offspring, calls] = evaluate(offspring, Task, 1);
            else
                calls = 0;
            end
        end
    end
end
