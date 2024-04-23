classdef riekf_SE3 < handle
    % Right-Invariant Extended Kalman filter class for 2D Localization, SE(2).
    %
    %   Author: Maani Ghaffari Jadidi
    %   Date:   01/23/2020
    
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
           R = X(1:3,1:3);
           RR = blkdiag(R,R);
           v = skew(X(1:3,4));
           p = skew(X(1:3,5));
           B = [v*R;p*R];
           %AdX = [X(1:3,1:3), [X(2,3); -X(1,3)]; 0 0 1];
           AdX = [R,zeros(3,6);
                 B ,RR];
        end
        
        function xhat = wedge(obj,x)
            % wedge operation for se(2) to put an R^3 vector into the Lie
            % algebra basis.
            Gox = [0     0     0   0   0;
                   0     0     -1   0   0;
                   0     1     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % omega_x
            
            Goy = [0     0     1   0   0;
                   0     0     0   0   0;
                   -1     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % omega_y

            Goz = [0     -1     0   0   0;
                   1     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % omega_z
%%%%%%
            Gvx = [0     0     0   1   0;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % v_x
            
            Gvy = [0     0     0   0   0;
                   0     0     0   1   0;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % v_y

            Gvz = [0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   1   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % v_z
%%%%%
            Gpx = [0     0     0   0   1;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % p_x

            Gpy = [0     0     0   0   0;
                   0     0     0   0   1;
                   0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   0]; % p_y

            Gpz = [0     0     0   0   0;
                   0     0     0   0   0;
                   0     0     0   0   1;
                   0     0     0   0   0;
                   0     0     0   0   0]; % p_z

            xhat = Gox * x(1) + Goy * x(2) + Goz * x(3) +  Gvx * x(4) + Gvy * x(5) + Gvz* x(6)+  Gpx * x(7) + Gpy * x(8) + Gpz * x(9);
        end
        
        function prediction(obj, u)
            % EKF propagation (prediction) step
            obj.P = obj.A * obj.P * obj.A' + ...
                obj.Ad(obj.X) * obj.Q * obj.Ad(obj.X)';
            % obj.P = obj.A * obj.P  + obj.P * obj.A' +...
            %    obj.Ad(obj.X) * obj.Q * obj.Ad(obj.X)';
          %dt = 0.001;
         % phi = expm(obj.A*dt);
         % Qr =  phi *obj.Ad(obj.X)*obj.Q*obj.Ad(obj.X)' * phi';

       %  obj.P = phi* obj.P *phi' + ...
        %       Qr*dt;
        obj.P
          obj.X = obj.f(obj.X,u);
        end
        
        function correction(obj, Y, b)
  % RI-EKF correction step
            %H = [obj.H(b1); obj.H(b2)]; % stack H
            %H = H([1:2,4:5],:); % 4x3 matrix, remove zero rows 
            H = obj.H;
            %N = obj.X * blkdiag(obj.N,0) * obj.X'; 
            %N = blkdiag(N(1:2,1:2), N(1:2,1:2)); % 4x4 block-diagonal matrix
            R = obj.X(1:3,1:3);

            s = 0.15;
            M = diag([0.15^(-2),0.15^(-2),0.15^(-2)]');

          %  M = diag(R*[0.15^2,0.15^2,0.15^2]');

             Sigma = inv(H*obj.Ad((obj.X))*obj.P*obj.Ad((obj.X))'*H');
             S_inv = Sigma - Sigma*inv(R'*M*R + Sigma)*Sigma;
             K = obj.P*(obj.Ad(obj.X)')*(H')*S_inv;

           % S = H * obj.P * H' + M;
           % K = (obj.P * H') * (S \ eye(size(S)))
            z = Y-obj.X(:,5);
           % z(1)=0;
            %z(2)=0;
            obsv = inv(obj.X)* (z);
            obsv = obsv(1:3);

            delta = obj.wedge(K*obsv);
            obj.X = expm(delta)*obj.X;
            % filter gain
            % S = H * obj.P * H' + N;
            % L = (obj.P * H') * (S \ eye(size(S)))

            I = eye(size(obj.P));
           % obj.P = (I - K * H) * obj.P * (I - K * H)' + K * M * K'; 
           % obj.P = (I - K * H) * obj.P * (I - K * H)' + K * M * K'; 
            %obj.P = 0.5*(obj.P+obj.P');
            obj.P = (eye(size(obj.P)) - K*H*obj.Ad((obj.X)))*obj.P;
            %b selects the position vector, while H selects just the z
            %measurement from it
           %  %Update State
           % nu = (blkdiag(obj.X, obj.X) * [Y1; Y2] - [b1; b2]);
           % Y(1) = obj.X(1,5);
           % Y(2) = obj.X(1,5);
           % Y(3) = Y(3);
           %  Yd = R'*Y(1:3);
           % 
           %  nu = obj.X*[Yd;0;1]; 
           % % nu(1)= 0;
           %  %nu(2) = 0;
           % 
           % % idx2keep_rows    = sum(abs(nu),2)>0 ;
           %  nu = nu(1:3)
           % % nu([3,6]) = [];
           %  delta = obj.wedge( L * nu) % innovation in the spatial frame
           %  obj.X = obj.f(obj.X,delta);
           % 
           %  % Update Covariance
           %  I = eye(size(obj.P));
           %  obj.P = (I - L * H) * obj.P * (I - L * H)' + L * N * L'; 


            % RI-EKF correction step
           %  %H = [obj.H(b1); obj.H(b2)]; % stack H
           %  %H = H([1:2,4:5],:); % 4x3 matrix, remove zero rows 
           %  H = obj.H;
           %  %N = obj.X * blkdiag(obj.N,0) * obj.X'; 
           %  %N = blkdiag(N(1:2,1:2), N(1:2,1:2)); % 4x4 block-diagonal matrix
           % 
           %  R = obj.X(1:3,1:3);
           %  N = diag(R*[0.15^2,0.15^2,0.15^2]');
           % % N = R*N*R';
           %  % filter gain
           %  S = H * obj.P * H' + N;
           %  L = (obj.P * H') * (S \ eye(size(S)))
           % 
           %  %%b selects the position vector, while H selects just the z
           %  %%measurement from it
           %  % Update State
           % % nu = (blkdiag(obj.X, obj.X) * [Y1; Y2] - [b1; b2]);
           % Y(1) = obj.X(1,5);
           % Y(2) = obj.X(1,5);
           % Y(3) = Y(3);
           %  Yd = R'*Y(1:3);
           % 
           %  nu = obj.X*[Yd;0;1]; 
           % % nu(1)= 0;
           %  %nu(2) = 0;
           % 
           % % idx2keep_rows    = sum(abs(nu),2)>0 ;
           %  nu = nu(1:3)
           % % nu([3,6]) = [];
           %  delta = obj.wedge( L * nu) % innovation in the spatial frame
           %  obj.X = obj.f(obj.X,delta);
           % 
           %  % Update Covariance
           %  I = eye(size(obj.P));
           %  obj.P = (I - L * H) * obj.P * (I - L * H)' + L * N * L'; 
        end



    end
end
