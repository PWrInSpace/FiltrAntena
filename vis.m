function  vis(X,Y,Z,omega,angle)

figure
subplot(2,1,1)
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
subplot(2,1,2)

hold on
grid on
plot(angle,omega)
legend("omega_x","omega_y","omega_z")
xlabel("angle[rad]");
ylabel("angular velocity[rad/s]");

end