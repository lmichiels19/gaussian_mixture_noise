
classdef NonLinearSystem < handle
    %SYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = NonLinearSystem()
            %SYSTEM Construct an instance of this class
        end
    end
    
    methods(Abstract)
        predict(obj,dT,x,u);
        measurement(obj,x);
    end
end

