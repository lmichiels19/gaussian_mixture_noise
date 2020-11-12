classdef LinearSystem < handle
    %SYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Cw_ % system noise
        Cv_ % measurement noise
    end
    
    methods
        function obj = LinearSystem(Cw,Cv)
            %SYSTEM Construct an instance of this class
            obj.Cw_ = Cw;
            obj.Cv_ = Cv;
        end
    end
    
    methods(Abstract)
        H(obj);
        A(obj);
        B(obj);
    end
    
    methods(Access=public)
        function Cw_ = Cw(obj)
            Cw_ = obj.Cw_;
        end
        
        function Cv_ = Cv(obj)
            Cv_ = obj.Cv_;
        end
    end
end
