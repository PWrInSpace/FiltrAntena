% Right-Invariant EKF localization on SE(2). The process model 
% is simply X_k+1 = X_k * exp(u_k + w_k) where
% X_k is in SE(2), u_k is the twist in se(2), and w_k is N(0,Q_k) and defined 
% in the Lie algebra se(2). The measurements are noisy 2D coordinates of the
% landmarks in Cartesian plane. We use expm and logm as numerical Lie exp 
% and log map. Both maps have closed-form formulas as well.
%
%   Author: Maani Ghaffari Jadidi
%   Date:   03/12/2020

clc; clear; close all

% incremental visualization
green = [0.2980 .6 0];
crimson = [220,20,60]/255; 
darkblue = [0 .2 .4];
Darkgrey = [.25 .25 .25];
VermillionRed = [156,31,46]/255;
DupontGray = [144,131,118]/255;

fsize = 14; % font size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');



X0T = [0; 0 ;5];%% initial position of the virtual target
X0M = [0 0 0];% initial position of the missile
T = 2*25;

V_M = 15; %%speed of the missile

V_T = 17; %% speed of the target

K = V_T/V_M; % only the ratio maters for determining the behaviour of the system
%K = 1.; %% if the ratio is < 1, the missile will catch up to the target, which is not what we want, so K >= 1
Kp = 1.5;
V0M = [6 1 9];% direction the missile's velocity

V0M = V0M/sqrt(V0M*V0M'); %% normalize since the ratio

out=sim('./trajectory_m.slx');

subplot(2,1,1)
grid on
hold on
plot3(out.Missile.Data(:,1),out.Missile.Data(:,2),out.Missile.Data(:,3))
l = length(out.Missile.Data(:,1));

sel = 1:floor(l/10):l;
quiver3(out.Missile.Data(sel,1),out.Missile.Data(sel,2),out.Missile.Data(sel,3), ...
    out.Vel.Data(sel,1),out.Vel.Data(sel,2),out.Vel.Data(sel,3),'-or')
xlabel('$x$', 'fontsize', fsize, 'Interpreter','latex')
ylabel('$y$', 'fontsize', fsize, 'Interpreter','latex')
zlabel('$z$', 'fontsize', fsize, 'Interpreter','latex')
%%
H = cell(length(out.Missile.Data(:,1)),1);
H{1} = eye(5);%Initial conditions

    v =   [out.Vel.Data(1,1),out.Vel.Data(1,2),out.Vel.Data(1,3)]';
    p =   [out.Missile.Data(1,1),out.Missile.Data(1,2),out.Missile.Data(1,3)]';
    v = v/vecnorm(v,2);
    roll = -atan2(v(2),v(3));
    pitch = asin(v(1));

    Rx = [1 0 0 ; 0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
    Ry = [cos(pitch) 0 sin(pitch); 0 1 0;-sin(pitch) 0 cos(pitch) ];
    
    R = Rx*Ry;
    dupa = posemat_SE2_3(R, v*V_M,p);
    H{1} = dupa;
%%

for i = 2:length(out.Missile.Data(:,1))
    v =   [out.Vel.Data(i,1),out.Vel.Data(i,2),out.Vel.Data(i,3)]';
    p =   [out.Missile.Data(i,1),out.Missile.Data(i,2),out.Missile.Data(i,3)]';

    roll = -atan2(v(2)/vecnorm(v,2),v(3)/vecnorm(v,2));
    pitch = asin(v(1)/vecnorm(v,2));

    Rx = [1 0 0 ; 0 cos(roll) -sin(roll);0 sin(roll) cos(roll)];
    Ry = [cos(pitch) 0 sin(pitch); 0 1 0;-sin(pitch) 0 cos(pitch) ];
    
    R = Rx*Ry;
    H{i} = posemat_SE2_3(R, v,p);
end
% generate noise-free twist control inputs (velocity commands) in the Lie algebra
u = cell(length(out.Missile.Data(:,1)),1);
u{1} = zeros(5);
for i = 2:length(out.Missile.Data(:,1))
    u{i} = logm(H{i-1} \H{i});
end


% construct noise free motion trajectory (sanity check for the generated
% inputs!)
path = [];
path.T = H{1};
path.x = 0;
path.y = 0;
path.z = 0;
path.vx = V0M(1)*V_M;
path.vy = V0M(2)*V_M;
path.vz = V0M(3)*V_M;
for i = 2:size(u,1)
   % path.T = path.T *expm(u{i})
    path.T = fx(path.T,u{i});
    path.x(i+1) = path.T(1,5);
    path.y(i+1) = path.T(2,5);
    path.z(i+1) = path.T(3,5);
    path.vx(i+1) = path.T(1,4);
    path.vy(i+1) = path.T(2,4);
    path.vz(i+1) = path.T(3,4);
end



subplot(2,1,2)
grid on
hold on
plot3(path.x,path.y,path.z)
sel = 1:floor(length(path.x)/10):length(path.x);
quiver3(path.x(sel),path.y(sel),path.z(sel), ...
        path.vx(sel),path.vy(sel),path.vz(sel), '-or')
xlabel('$x$', 'fontsize', fsize, 'Interpreter','latex')
ylabel('$y$', 'fontsize', fsize, 'Interpreter','latex')
zlabel('$z$', 'fontsize', fsize, 'Interpreter','latex')

 % build a system 
 sys = [];
 % motion model noise covariance
 sys.X = dupa;
 sys.Q = diag([repmat(0.015^2,1,3), repmat(0.005^2,1,3), repmat(0,1,3)]);
 %sys.A = eye(9);
 sys.A = [zeros(6,9);
          zeros(3,3), eye(3),zeros(3,3)];
 %sys.f = @(x,u) x*expm(u);
sys.f = @(x,u) fx(x,u);


 %sys.H = @(m) [m(2) -1 0; -m(1) 0 -1; 0 0 0]; % this is a cuntion to serve
 %as a siplification in constructing the actual H matrix in the filter
 %class, a constant H will suffice for barometer measurements
 sys.H = [zeros(3,6),eye(3)];
 sys.N = 0.15^2;%We measure 1 value, so the matrix is just a constant
% sys.N = diag([0.5^2; 0.5^2]);
% 
% % se(2) Lie algebra basis twist = vec(\omega, v_1, v_2)
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
                   0     0     0   0   0
                   ]; % p_z
% % now make the twist noisy! in practice the velocity readings are not
% % perfect.
% % Cholesky factor of covariance for sampling
%LQ = chol(sys.Q, 'lower');
u_n=cell(length(out.Missile.Data(:,1)),1);
u_n{1}=zeros(5);
for i = 1:length(u)
    noise = sys.Q* randn(9,1);
    %noise = [0,0,0];
    N = Gox * noise(1) + Goy * noise(2) + Goz * noise(3) + Gvx * noise(4) + Gvy * noise(5) + Gvz * noise(6) + Gpx * noise(7) + Gpy * noise(8) + Gpz * noise(9);
    u_n{i} = u{i} + N;
end

path_noisy = [];
path_noisy.T = H{1};
path_noisy.x = 0;
path_noisy.y = 0;
path_noisy.z = 0;
path_noisy.vx = 0;
path_noisy.vy = 0;
path_noisy.vz =0;
acc=zeros(3,l);
acc(:,1) = [0;0;0];
for i = 1:size(u,1)
    %path_noisy.T =   path_noisy.T *expm(u_n{i});
    path_noisy.T =   fx(path_noisy.T,u_n{i});
    path_noisy.x(i+1) = path_noisy.T(1,5);
    path_noisy.y(i+1) = path_noisy.T(2,5);
    path_noisy.z(i+1) = path_noisy.T(3,5);
    path_noisy.vx(i+1) = path_noisy.T(1,4);
    path_noisy.vy(i+1) = path_noisy.T(2,4);
    path_noisy.vz(i+1) = path_noisy.T(3,4);
    %% testing which reference frame this applies to
     %acc(:,i) = u{i}(1:3,4);
     % acc(3,i) = -u{i}(1,2);
     % acc(2,i) = u{i}(1,3);
     % acc(1,i) = -u{i}(2,3);
end


% figure 
% subplot(2,1,1)
% grid on
% hold on
% 
% plot(1:l,out.Acc.Data)
% 
% subplot(2,1,2)
% grid on
% hold on
% 
% plot(1:l,acc)



subplot(2,1,2)
grid on
hold on
plot3(path_noisy.x,path_noisy.y,path_noisy.z)
sel = 1:floor(length(path.x)/10):length(path.x);
quiver3(path_noisy.x(sel),path_noisy.y(sel),path_noisy.z(sel), ...
        path_noisy.vx(sel),path_noisy.vy(sel),path_noisy.vz(sel), '-og')
xlabel('$x$', 'fontsize', fsize, 'Interpreter','latex')
ylabel('$y$', 'fontsize', fsize, 'Interpreter','latex')
zlabel('$z$', 'fontsize', fsize, 'Interpreter','latex')

legend("True trajectory","True Velocity","Noise integrated trajecotry","Noise integrated velocity")


h_leg = []; % legend handle
figure; hold on
h_leg{1} = plot3(path.x, path.y,path.z,'-', 'color', Darkgrey, 'linewidth', 3);
grid on%, axis auto equal, axis([-6 6 -1 5])
xlabel('$x$', 'fontsize', fsize, 'Interpreter','latex')
ylabel('$y$', 'fontsize', fsize, 'Interpreter','latex')
zlabel('$z$', 'fontsize', fsize, 'Interpreter','latex')
set(gca, 'fontsize', fsize)
%axis equal, axis([-l/4 l*1.35 -l/4 l*1.35]), axis off

 filter = riekf_SE3(sys); % create an RI-EKF object

 % plot initial mean 
 h_leg{2} = plot3(filter.X(1,5), filter.X(2,5),filter.X(3,5), 'o', 'color', [crimson, 0.7], 'MarkerFaceColor', crimson, 'markersize', 8);
 h_leg{3} = quiver3(filter.X(1,5), filter.X(2,5),filter.X(3,5), ...
            filter.X(1,4), filter.X(2,4),filter.X(3,4),'Color','red');

% video recorder object
% video = VideoWriter('dead-reckoning_se2.avi', 'Motion JPEG AVI');
% view(90,0)
% 
% video = VideoWriter('riekf_localization_se2.avi', 'Motion JPEG AVI');
% video.Quality = 100;
% open(video)

skip = 10000;

path_filter = [];
path_filter.T = H{1};
path_filter.x = 0;
path_filter.y = 0;
path_filter.z = 0;
path_filter.vx = 0;
path_filter.vy = 0;
path_filter.vz =0;

 h_leg{4} = plot3(path_filter.x,path_filter.y,path_filter.z);

for i = 1:size(u,1)
    % predict next pose using given twist
    filter.prediction(u_n{i});

    if ~mod(i,skip)

        b = [0 0 0 0 1]';
        Y = H{i}*b ;%+ [ randn(3,1)/100; 0;0];

       % correction based on the measurements
        filter.correction(Y, b);
    end


    path_filter.x(i+1) = filter.X(1,5);
    path_filter.y(i+1) = filter.X(2,5);
    path_filter.z(i+1) = filter.X(3,5);
    path_filter.vx(i+1) = filter.X(1,4);
    path_filter.vy(i+1) = filter.X(2,4);
    path_filter.vz(i+1) = filter.X(3,4);

   % update graphics
    set(h_leg{2},'XData',filter.X(1,5),'YData', filter.X(2,5),'ZData', filter.X(3,5));
    set(h_leg{3},'XData',filter.X(1,5),'YData', filter.X(2,5),'ZData', filter.X(3,5), 'UData', filter.X(1,4), 'VData', filter.X(2,4), 'WData', filter.X(3,4));
    set(h_leg{4},'XData',path_filter.x,'YData', path_filter.y,'ZData', path_filter.z);

     %drawnow 
     %pause(0.0001)
    %  frame = getframe(gcf);
    % writeVideo(video,frame);
end


% figure
% hold on
% grid on
% plot3(path_filter.x,path_filter.y,path_filter.z)

% close(video)

figure
subplot(2,2,1)
hold on
grid on
plot(1:length(path.x),path.x-path_filter.x)
plot(1:length(path.x),path.y-path_filter.y)
plot(1:length(path.x),path.z-path_filter.z)
title("pos err filter")
subplot(2,2,2)
hold on
grid on
plot(1:length(path.x),path.x-path_noisy.x)
plot(1:length(path.x),path.y-path_noisy.y)
plot(1:length(path.x),path.z-path_noisy.z)
% figuresize(21,21,'cm')
% print -opengl -dpng -r600 riekf_loc_se2.png
title("pos err noise")

subplot(2,2,3)
hold on
grid on
plot(1:length(path.vx),path.vx-path_filter.vx)
plot(1:length(path.vy),path.vy-path_filter.vy)
plot(1:length(path.vz),path.vz-path_filter.vz)
title("vel err filter")
subplot(2,2,4)
hold on
grid on
plot(1:length(path.vx),path.vx-path_noisy.vx)
plot(1:length(path.vy),path.vy-path_noisy.vy)
plot(1:length(path.vz),path.vz-path_noisy.vz)
% figuresize(21,21,'cm')
% print -opengl -dpng -r600 riekf_loc_se2.png
title("vel err noise")

function H = posemat_SE_2(x,y,h)
% construct a SE(2) matrix element
H = [cos(h) -sin(h) x;
    sin(h)  cos(h)  y;
    0       0       1];
end



function ELLIPSE = confidence_ellipse(X,L)
% create confidence ellipse
% se(2) Lie algebra basis twist = vec(\omega, v_1, v_2)
G1 = [0    -1     0
      1     0     0
      0     0     0];
G2 = [0     0     1
      0     0     0
      0     0     0];
G3 = [0     0     0
      0     0     1
      0     0     0];

% first create points from a unit circle + angle (third dimension of so(3))
phi = (-pi:.01:pi)';
circle = [zeros(length(phi),1), cos(phi), sin(phi)];
% Chi-squared 3-DOF 95% percent confidence (0.05): 7.815
scale = sqrt(7.815);
% main loop; iterate over the control inputs and move the robot
ELLIPSE = zeros(size(circle,1),2); % covariance ellipse on manifold (nonlinear)

for j = 1:size(circle,1)
    % sample covariance on SE(2)
    ell_se2_vec = scale * L * circle(j,:)';
    % retract and left-translate the ellipse on Lie algebra to SE(2) using Lie exp map
    temp = X * expm(G1 * ell_se2_vec(1) + G2 * ell_se2_vec(2) + G3 * ell_se2_vec(3));
    ELLIPSE(j,:) = [temp(1,3), temp(2,3)];
end



end

function fx = fx(x,u)
    f = expm(u);
    f(1:3,5) = [0 0 0]';
    R = x(1:3,1:3);
    Z = zeros(5);
    dt = 0.01;
    Z(1:3,5) =  x(1:3,4)*dt+R*f(1:3,4)*dt/2;
    fx = x*f + Z;
end