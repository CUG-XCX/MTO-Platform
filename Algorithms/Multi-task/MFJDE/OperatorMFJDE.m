classdef OperatorMFJDE < OperatorDE
    methods (Static)
        function [offspring, calls] = generate(callfun, population, Tasks, rmp, t1, t2)
             if length(population) <= 3
                offspring = population;
                calls = 0;
                return;
            end
            Individual_class = class(population(1));
            group = cell([1, length(Tasks)]);
            for i = 1:length(population)
                group{population(i).skill_factor} = [group{population(i).skill_factor}, i];
            end
           for i = 1:length(population)
                offspring(i) = feval(Individual_class);
                offspring(i).factorial_costs = inf(1, length(Tasks));
                other = [];
                for t = 1:length(group)
                    if population(i).skill_factor ~= t
                        other = [other, group{t}];
                    end
                end
                other = other(randperm(length(other)));
                A = randperm(length(group{population(i).skill_factor}));
                A = group{population(i).skill_factor}(A);
                A(A == i) = [];
                N = length(population);
                idx = 1;
                
                if rand < rmp
                    x1 = A(1);
                    x2 = other(mod(2 - 1, length(other)) + 1);
                    x3 = other(mod(3 - 1, length(other)) + 1);
                    offspring(i).skill_factor = population(x2).skill_factor;
                else
                    %rank_jde
                   while rand > (N - population(A(idx)).factorial_ranks(population(i).skill_factor))/ N || A(idx) == i
                     idx = mod(idx, length(A)) + 1;
                   end
                   x1 = A(idx);
                   while rand > (N - population(A(idx)).factorial_ranks(population(i).skill_factor)) / N || A(idx) == x1 || A(idx) == i
                     idx = mod(idx, length(A)) + 1;
                   end
                   x2 = A(idx);
                   while A(idx) == x1 || A(idx) == x2 || A(idx) == i
                     idx = mod(idx, length(A)) + 1;
                   end
                   x3 = A(idx); 
                   offspring(i).skill_factor = population(x1).skill_factor;
                end   
              
        
                % parameter self-adaptation
                offspring(i).F = population(i).F;
                offspring(i).CR = population(i).CR;
                if rand < t1
                    offspring(i).F = rand * 0.9 + 0.1;
                end
                if rand < t2
                    offspring(i).CR = rand;
                end

                offspring(i) = OperatorjDE_rank.mutate_rand_1(offspring(i), population(x1), population(x2), population(x3));
                offspring(i) = OperatorjDE_rank.crossover(offspring(i), population(i));

                offspring(i).rnvec(offspring(i).rnvec > 1) = 1;
                offspring(i).rnvec(offspring(i).rnvec < 0) = 0;
                %disp("x");
            end
            if callfun
                offspring_temp = feval(Individual_class).empty();
                calls = 0;
                for t = 1:length(Tasks)
                    offspring_t = offspring([offspring.skill_factor] == t);
                    [offspring_t, cal] = evaluate(offspring_t, Tasks(t), t);
                    offspring_temp = [offspring_temp, offspring_t];
                    calls = calls + cal;
                end
                offspring = offspring_temp;
            else
                calls = 0;
            end
           % disp("1");
        end
    end
end