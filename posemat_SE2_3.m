
function H = posemat_SE2_3(R,v,x)
% construct a SE(2) matrix element
H = [R,v,x;
     zeros(2,3),eye(2)];
end