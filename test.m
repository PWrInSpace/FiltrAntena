clear all
close all
clc




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

function mat = Rod(theta,x)
skew = [0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
mat = eye(3) + sin(theta)*skew + (1-cos(theta))*skew^2;
end