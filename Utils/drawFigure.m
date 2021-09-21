classdef drawFigure < handle

    properties
        x_cell
        y_cell
        xlabel_str
        ylabel_str
        legend_cell
        title_str
        save_dir = './'
        figure_type = 'png'
    end

    methods

        function obj = setXY(obj, x_cell, y_cell)
            obj.x_cell = x_cell;
            obj.y_cell = y_cell;
        end

        function obj = setXYlabel(obj, xlabel_str, ylabel_str)
            obj.xlabel_str = xlabel_str;
            obj.ylabel_str = ylabel_str;
        end

        function obj = setLegend(obj, legend_cell)
            obj.legend_cell = legend_cell;
        end

        function obj = setTitle(obj, title_str)
            obj.title_str = title_str;
        end

        function obj = setSaveDir(obj, save_dir)
            obj.save_dir = save_dir;
        end

        function obj = setFigureType(obj, figure_type)
            obj.figure_type = figure_type;
        end

        function obj = draw(obj)
            fig = figure('Visible', 'off');

            for i = 1:length(obj.x_cell)
                plot(obj.x_cell{i}, obj.y_cell{i}, ':', 'LineWidth', 3);
                hold on;
            end

            title(strrep(obj.title_str, '_', '\_'));
            xlabel(strrep(obj.xlabel_str, '_', '\_'));
            ylabel(strrep(obj.ylabel_str, '_', '\_'));
            legend(strrep(obj.legend_cell, '_', '\_'));
            max_x = 0;

            for i = 1:length(obj.x_cell)

                if obj.x_cell{i}(end) > max_x
                    max_x = obj.x_cell{i}(end);
                end

            end

            xlim([0, max_x]);

            file_name = [obj.save_dir, obj.title_str, '.', obj.figure_type];

            if strcmp(obj.figure_type, 'eps')
                saveas(fig, file_name, 'psc2');
            else
                saveas(fig, file_name);
            end

        end

    end

end
