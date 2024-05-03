function mat = Rod(theta,x)
skew = [0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
mat = eye(3) + sin(theta)*skew + (1-cos(theta))*skew^2;
end