function [x_grid] = create_sample_grid(limits, count)
%CREATE_SAMPLE_GRID Grid of n dimensions with limits and point count
% 1 to 4 Dimensions are supported

    % Dimension count
    D = size(limits,1);
    
    if D == 1
        x_grid = linspace(limits(1,1),limits(1,2),count(1));
    elseif D == 2
        d1 = linspace(limits(1,1),limits(1,2),count(1));
        d2 = linspace(limits(2,1),limits(2,2),count(2));
        [D1, D2] = ndgrid(d1,d2);
        x_grid = [D1(:),D2(:)]';
    elseif D == 3
        d1 = linspace(limits(1,1),limits(1,2),count(1));
        d2 = linspace(limits(2,1),limits(2,2),count(2));
        d3 = linspace(limits(3,1),limits(3,2),count(3));
        [D1, D2, D3] = ndgrid(d1,d2,d3);
        x_grid = [D1(:),D2(:),D3(:)]';
    elseif D == 4
        d1 = linspace(limits(1,1),limits(1,2),count(1));
        d2 = linspace(limits(2,1),limits(2,2),count(2));
        d3 = linspace(limits(3,1),limits(3,2),count(3));
        d4 = linspace(limits(4,1),limits(4,2),count(4));
        [D1, D2, D3, D4] = ndgrid(d1,d2,d3,d4);
        x_grid = [D1(:),D2(:),D3(:),D4(:)]';
    else
        nd_grid_parameter = "";
        nd_grid_output = "[";
        x_gird_cm = "x_grid = [";
        for i=1:D
            eval("d"+i+" = linspace(limits("+i+",1),limits("+i+",2),count("+i+"));");
            if i ~= D
                nd_grid_parameter = nd_grid_parameter+"d"+i+",";
                nd_grid_output = nd_grid_output+"D"+i+",";
                x_gird_cm = x_gird_cm+"D"+i+"(:),";
            else   % no comma seperator
                nd_grid_parameter = nd_grid_parameter+"d"+i; 
                nd_grid_output = nd_grid_output+"D"+i;
                x_gird_cm = x_gird_cm+"D"+i+"(:)";
            end
        end
        x_gird_cm = x_gird_cm+"]';";
        nd_grid_output = nd_grid_output+"]";
        eval(nd_grid_output+" = ndgrid("+nd_grid_parameter+");");
        eval(x_gird_cm);
    end

end

