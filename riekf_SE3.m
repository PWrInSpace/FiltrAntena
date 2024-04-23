classdef riekf_SE3 < handle
    % Right-Invariant Extended Kalman filter class for orientation estumation SO(3).
    %
    %   Authors: Jacek Grzegorzewski, Patryk Niczke
    %   Date:  23.04.204
    %   Based on work by: Maani Ghaffari Jadidi
    %   Link to original repo: to be added
    
    properties
        A;              % error dynamics matrix
        H;              % measurement error matrix
        f;              % process model
        X;              % state vector
        P;              % state covariance
        Q;              % input noise covariance
        N;              % measurement noise covariance
        NU;
    end
    
    methods
        function obj = riekf_SE3(system)
            % riekf Construct an instance of this class
            %
            %   Inputs:
            %       system          - system and noise models
            obj.A = system.A;
            obj.f = system.f;
            obj.H = system.H;
            obj.Q = system.Q;
            obj.N = system.N;
            %obj.X = eye(5);% state matrix is 5x5
            obj.X = system.X;
            obj.P = 0.5*ones(9); %% P's dimension depends on te dimension of the lie algebra
        end
        
        function AdX = Ad(obj, X)
           % Adjoint
           %Adjoint of R in SO(3) is just R
           AdX = X;
        end
        
        function xhat = wedge(obj,x)
            % wedge operation for se(2) to put an R^3 vector into the Lie
            % algebra basis.
            Gox = [0     0     0  ;
                   0     0     -1  ;
                   0     1     0   ]; % omega_x
            
            Goy = [0     0     1   ;
                   0     0     0   ;
                   -1     0     0  ]; % omega_y

            Goz = [0     -1     0  ;
                   1     0     0   ;
                   0     0     0   ]; % omega_z
            xhat = Gox * x(1) + Goy * x(2) + Goz * x(3) ;
        end
        
        function prediction(obj, u)
            %TODO
            %Implement right invariant prediction with constant A matrix
        end
        
        function correction_IMU(obj, Y, b)
            %TODO
            %Implement accelereometer based correction step
        end

        function correction_MAG(obj, Y, b)
            %TODO
            %Implement magnetometer based correction step
        end


    end
end
