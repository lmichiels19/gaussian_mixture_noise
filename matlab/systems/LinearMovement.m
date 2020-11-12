classdef LinearMovement < NonLinearSystem & LinearSystem
    %BALLSYSTEM Non-linear System for a Ball flying in a circle
    % System state (x)
    % x - x position (m)
    % y - y position (m)
    % Measurement (z)
    % x - position
    % y - position
    % Control (u)
    % x - speed
    % y - speed
    
    properties
    end
    
    methods
        function obj = LinearMovement(Cw,Cv)
            obj = obj@LinearSystem(Cw,Cv);
            %SYSTEM Construct an instance of this class
        end
                
        
        function x = predict(obj,x,dT,u)
            %predict Predicts system state
            % \param dT - timestep
            % \param u - command vector

            if size(x,2) > 1
                % list of states
                x(1,:) = x(1,:)+x(3,:).*dT;
                x(2,:) = x(2,:)+x(4,:).*dT;
                x(3,:) = x(3,:)+u(1);
                x(4,:) = x(4,:)+u(2);
            else           
                x(1) = x(1)+x(3).*dT;
                x(2) = x(2)+x(4).*dT;
                x(3) = x(3)+u(1);
                x(4) = x(4)+u(2);
            end
        end
        
        function z = measurement(obj,x)
            %measurement Returns measurement for system state
            z = zeros(2,1);
            z = x;
        end
        
        function z = measurementSample(obj,x)
            %measurement Returns measurement for system state
            z = mvnrnd(x(1:2)',obj.Cv)';
        end
        
        function w = systemNoiseSample(obj,x,dT)
            if size(x,2) > 1
                % list of noise samples
                w = mvnrnd(zeros(size(x,1),1),obj.Cw,size(x,2))';
            else
                % on noise sample
                w = mvnrnd(zeros(size(x,1),1),obj.Cw)';
            end
        end
        
        function li = likelihood(obj,x,z)
            N = size(x,2);
            if N > 1
                li = mvnpdf(z',x(1:2,:)',obj.Cv)'.*sqrt(2*pi*det(obj.Cv));
            else
                li = mvnpdf(z',x(1:2)',obj.Cv)'.*sqrt(2*pi*det(obj.Cv));
            end
        end
        
        function A_ = A(obj,dT)
            A_ = [1 0 dT 0; 0 1 0 dT; 0 0 1 0; 0 0 0 1];
        end
        
        function B_ = B(obj,dT)
            B_ = [0 0 dT dT];
        end
        
        function H_ = H(obj)
            H_ = [1 0 0 0; 0 1 0 0];
        end
    end
end

