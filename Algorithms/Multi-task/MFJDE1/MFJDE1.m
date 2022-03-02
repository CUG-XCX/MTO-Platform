classdef MFJDE1 < Algorithm
    % <Multi> <None>
    properties (SetAccess = private)
        rmp = 0.3
        t1 = 0.1;
        t2 = 0.1;
    end

    methods
        function parameter = getParameter(obj)
            parameter = {'rmp: Random Mating Probability', num2str(obj.rmp), ...
                        't1: probability of F change', num2str(obj.t1), ...
                        't2: probability of CR change', num2str(obj.t2)};
        end

        function obj = setParameter(obj, parameter_cell)
            count = 1;
            obj.rmp = str2double(parameter_cell{count}); count = count + 1;
            obj.t1 = str2double(parameter_cell{count}); count = count + 1;
            obj.t2 = str2double(parameter_cell{count}); count = count + 1;
        end

        function data = run(obj, Tasks, run_parameter_list)
            sub_pop = run_parameter_list(1);
            iter_num = run_parameter_list(2);
            eva_num = run_parameter_list(3) * length(Tasks);
            tic

            pop_size = sub_pop * length(Tasks);

            % initialize
            [population, fnceval_calls, bestobj, data.bestX] = initializeMF(IndividualjDE, pop_size, Tasks, length(Tasks));
    
            for i = 1:length(population)  
                % initialize F and CR
                 population(i).F = rand * 0.9 + 0.1;
                 population(i).CR = rand;
                 %disp(population(i).F);
                 %disp(population(i).factorial_ranks);
            end
            data.convergence(:, 1) = bestobj;
            generation = 1;
            while generation < iter_num && fnceval_calls < eva_num
                generation = generation + 1;
                
                % generation
                [offspring, calls] = OperatorMFJDE.generate(1, population, Tasks, obj.rmp, obj.t1, obj.t2);
                fnceval_calls = fnceval_calls + calls;

                % selection
                [population, bestobj, data.bestX] = selectMF(population, offspring, Tasks, pop_size, bestobj, data.bestX);
                data.convergence(:, generation) = bestobj;
               
                
            end
            data.bestX = uni2real(data.bestX, Tasks);
            %disp("3");
            data.clock_time = toc;
        end
    end
end