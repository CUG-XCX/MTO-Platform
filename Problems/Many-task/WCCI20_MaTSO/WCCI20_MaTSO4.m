classdef WCCI20_MaTSO4 < Problem
    % <Many> <None>

    properties
        task_size = 50;
        dims = 50;
    end

    methods
        function parameter = getParameter(obj)
            parameter = {'Task Num', num2str(obj.task_size), ...
                        'Dims', num2str(obj.dims)};
            parameter = [obj.getRunParameter(), parameter];
        end

        function obj = setParameter(obj, parameter_cell)
            obj.setRunParameter(parameter_cell(1:2));
            count = 3;
            obj.task_size = str2double(parameter_cell{count}); count = count + 1;
            obj.dims = str2double(parameter_cell{count}); count = count + 1;
        end

        function Tasks = getTasks(obj)
            Tasks = benchmark_WCCI20_MaTSO(4, obj.task_size, obj.dims);
        end

    end

end
