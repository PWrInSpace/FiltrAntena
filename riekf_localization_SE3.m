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

angle = 0:0.001:10;

v = @(x) [cos(x),sin(x),0];
H = cell(length(angle));
u = cell(length(angle));
H{1} = eye(3);%Initial conditions

X = [ones(1,length(angle));
     zeros(2,length(angle))];

Y = [zeros(1,length(angle));
     ones(1,length(angle));
     zeros(1,length(angle))];

Z = [zeros(2,length(angle));
    ones(1,length(angle))];

omega = zeros(3,length(angle));

for i = 1:length(angle)

    H{i} = Rod((angle(i))*2*pi*0.5  ,v(angle(i)*2*pi ));
    X(:,i) = H{i}*X(:,i);
    Y(:,i) = H{i}*Y(:,i);
    Z(:,i) = H{i}*Z(:,i);
    if(i > 1)

        u{i} = logm(H{i-1}\H{i});
        omega(:,i) = [-u{i}(2,3),u{i}(1,3), -u{i}(1,2)];
    end
end


figure
grid on
hold on
axis([-1 1 -1 1 -1 1])
plot3(X(1,:),X(2,:),X(3,:),'r')
plot3(Y(1,:),Y(2,:),Y(3,:),'g')
plot3(Z(1,:),Z(2,:),Z(3,:),'b')

 quiver3(0,0,0, ...
           1,0,0,'Color','red');
 quiver3(0,0,0, ...
           0,1,0,'Color','green');
  quiver3(0,0,0, ...
           0,0,1,'Color','blue');

  xlabel("X")
  ylabel("Y")
  zlabel("Z")
legend("X","Y","Z")

figure
hold on
grid on
plot(angle,omega)


len = length(angle)



H = cell(len);
%%TODO add initial conditions
H{1} = eye(3);
%%
%%Generate ground truth orientations
for i = 2:length(len)
    H{i} = ...
end

% generate noise-free twist control inputs (velocity commands) in the Lie algebra
u = cell(len)
u{1} = zeros(3);
for i = 1:len
    u{i} = logm(H{i-1} \H{i});
end



% construct noise free motion trajectory (sanity check for the generated
% inputs!)
path = [];
path.T = H{1};
for i = 2:size(u,1)
    path.T = path.T *expm(u{i})
end


 %TODO
 % build a system 
 sys = [];
 % motion model noise covariance
 sys.X = ...;
 sys.Q = ...;
 sys.A = ...;
 sys.f = ...;
 sys.H = ...;
 sys.N = ...;
% % so(3) Lie algebra basis twist 
            Gox = [0     0     0  ;
                   0     0     -1  ;
                   0     1     0   ]; % omega_x
            
            Goy = [0     0     1   ;
                   0     0     0   ;
                   -1     0     0  ]; % omega_y

            Goz = [0     -1     0  ;
                   1     0     0   ;
                   0     0     0   ]; % omega_z
 

% % now make the twist noisy! in practice the velocity readings are not
% % perfect.
u_n=cell(len);
u_n{1}=zeros(3);
for i = 1:length(u)
    noise = sys.Q* randn(9,1);
    %noise = [0,0,0];
    N = Gox * noise(1) + Goy * noise(2) + Goz * noise(3);
    u_n{i} = u{i} + N;
end

path_noisy = [];
path_noisy.T = H{1};
for i = 1:size(u,1)
    %TODO 
    %integrate noisy measurements
    path_noisy.T =  ...;
end




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

%TODO 
%choose measurement frequencies
skip_mag = 10000;
skip_acc = 10000;

path_filter = [];
path_filter.T = H{1};

 h_leg{4} = plot3(path_filter.x,path_filter.y,path_filter.z);

for i = 1:size(u,1)
    % predict next pose using given twist
    filter.prediction(u_n{i});

    if ~mod(i,skip_mag)
        %TODO 
        %implement magnetometer correction step
    end

    if ~mod(i,skip_acc)
        %TODO
        %implement accelerometer correction step
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

% close(video)

% figuresize(21,21,'cm')
% print -opengl -dpng -r600 riekf_loc_se2.png
% figuresize(21,21,'cm')
% print -opengl -dpng -r600 riekf_loc_se2.png

function H = posemat_SE_2(x,y,h)
% construct a SE(2) matrix element
H = [cos(h) -sin(h) x;
    sin(h)  cos(h)  y;
    0       0       1];
end




function fx = fx(x,u)
    %%TODO
    %define state change function
    %e.g. system dynamics, x_{k+1} = f(x_k,u)
end


