syms q1 q2 q3 q4 q5 q6 x1 y1 z1 x2 y2 z2 f u_0 v_0 u01 v01 u02 v02;
Rot_x=[1 0 0; 0 cos(q1) sin(q1); 0 -sin(q1) cos(q1)];   % Rotation matrix over x axis
Rot_y=[cos(q2) 0 -sin(q2); 0 1 0; sin(q2) 0 cos(q2)];   % Rotation matrix over y axis
Rot_z=[cos(q3) sin(q3) 0; -sin(q3) cos(q3) 0; 0 0 1];   % Rotation matrix over z axis
rot_map=Rot_z*Rot_y*Rot_x;      % Eauler rotation map
trans_vec=[q4; q5; q6];         % Translation vector of projection matrix
Proj_mat1=[rot_map trans_vec]*[x1; y1; z1;1];
Proj_mat2=[rot_map trans_vec]*[x2; y2; z2;1];
u_coordinate1=-f*Proj_mat1(2,1)+(u_0-u01)*Proj_mat1(3,1);
v_coordinate1=f*Proj_mat1(1,1)+(v_0-v01)*Proj_mat1(3,1);
u_coordinate2=-f*Proj_mat2(2,1)+(u_0-u02)*Proj_mat2(3,1);
v_coordinate2=f*Proj_mat2(1,1)+(v_0-v02)*Proj_mat2(3,1);

u_coordinate1=Proj_mat1(2,1)-(u01)*Proj_mat1(3,1);
v_coordinate1=Proj_mat1(1,1)-(v01)*Proj_mat1(3,1);
u_coordinate2=Proj_mat2(2,1)-(u02)*Proj_mat2(3,1);
v_coordinate2=Proj_mat2(1,1)-(v02)*Proj_mat2(3,1);

[-f (u_0-u01) -f*rot_map(2,:)+(u_0-u01)*rot_map(3,:)]
[ f (v_0-v01)  f*rot_map(1,:)+(v_0-v01)*rot_map(3,:)]
[-f (u_0-u02) -f*rot_map(2,:)+(u_0-u02)*rot_map(3,:)]
[ f (v_0-v02)  f*rot_map(1,:)+(v_0-v02)*rot_map(3,:)]

u1=u_0-f*Proj_mat1(2,1)/Proj_mat1(3,1);
v1=v_0+f*Proj_mat1(1,1)/Proj_mat1(3,1);
u2=u_0-f*Proj_mat2(2,1)/Proj_mat2(3,1);
v2=v_0+f*Proj_mat2(1,1)/Proj_mat2(3,1);

E=(u_coordinate1)^2+(v_coordinate1)^2+(u_coordinate2)^2+(v_coordinate2)^2;

G(1,1)=diff(E,q4);
G(1,2)=diff(E,q5);
G(1,3)=diff(E,q6);
G(1,4)=diff(E,x1);
G(1,5)=diff(E,y1);
G(1,6)=diff(E,z1);
G(1,7)=diff(E,x2);
G(1,8)=diff(E,y2);
G(1,9)=diff(E,z2);


H(1,1:6)=[diff(diff(E,x1),x1) diff(diff(E,x1),y1) diff(diff(E,x1),z1) diff(diff(E,x1),x2) diff(diff(E,x1),y2) diff(diff(E,x1),z2)];
H(2,2:6)=[diff(diff(E,y1),y1) diff(diff(E,y1),z1) diff(diff(E,y1),x2) diff(diff(E,y1),y2) diff(diff(E,y1),z2)];
H(3,3:6)=[diff(diff(E,z1),z1) diff(diff(E,z1),x2) diff(diff(E,z1),y2) diff(diff(E,z1),z2)];
H(4,4:9)=[diff(diff(E,x2),x2) diff(diff(E,x2),y2) diff(diff(E,x2),z2)];
H(5,5:9)=[diff(diff(E,y2),y2) diff(diff(E,y2),z2)];
H(6,6)=diff(diff(E,z2),z2);

G1(1,1)=diff(u_coordinate1,q1);
G1(1,2)=diff(u_coordinate1,q2);
G1(1,3)=diff(u_coordinate1,q3);
G1(1,4)=diff(u_coordinate1,q4);
G1(1,5)=diff(u_coordinate1,q5);
G1(1,6)=diff(u_coordinate1,q6);

G1(1,7)=diff(v_coordinate1,q1);
G1(1,8)=diff(v_coordinate1,q2);
G1(1,9)=diff(v_coordinate1,q3);
G1(1,10)=diff(v_coordinate1,q4);
G1(1,11)=diff(v_coordinate1,q5);
G1(1,12)=diff(v_coordinate1,q6);

G1(1,13)=diff(u_coordinate2,q1);
G1(1,14)=diff(u_coordinate2,q2);
G1(1,15)=diff(u_coordinate2,q3);
G1(1,16)=diff(u_coordinate2,q4);
G1(1,17)=diff(u_coordinate2,q5);
G1(1,18)=diff(u_coordinate2,q6);

G1(1,19)=diff(u_coordinate2,q1);
G1(1,20)=diff(u_coordinate2,q2);
G1(1,21)=diff(u_coordinate2,q3);
G1(1,22)=diff(u_coordinate2,q4);
G1(1,23)=diff(u_coordinate2,q5);
G1(1,24)=diff(u_coordinate2,q6);


J(1,1)=diff(u_coordinate1,q1);
J(1,2)=diff(u_coordinate1,q2);
J(1,3)=diff(u_coordinate1,q3);
J(1,4)=diff(u_coordinate1,q4);
J(1,5)=diff(u_coordinate1,q5);
J(1,6)=diff(u_coordinate1,q6);
J(1,7)=diff(u_coordinate1,x1);
J(1,8)=diff(u_coordinate1,y1);
J(1,9)=diff(u_coordinate1,z1);
J(1,10)=diff(u_coordinate1,x2);
J(1,11)=diff(u_coordinate1,y2);
J(1,12)=diff(u_coordinate1,z2);

J(2,1)=diff(v_coordinate1,q1);
J(2,2)=diff(v_coordinate1,q2);
J(2,3)=diff(v_coordinate1,q3);
J(2,4)=diff(v_coordinate1,q4);
J(2,5)=diff(v_coordinate1,q5);
J(2,6)=diff(v_coordinate1,q6);
J(2,7)=diff(v_coordinate1,x1);
J(2,8)=diff(v_coordinate1,y1);
J(2,9)=diff(v_coordinate1,z1);
J(2,10)=diff(v_coordinate1,x2);
J(2,11)=diff(v_coordinate1,y2);
J(2,12)=diff(v_coordinate1,z2);

J(3,1)=diff(u_coordinate2,q1);
J(3,2)=diff(u_coordinate2,q2);
J(3,3)=diff(u_coordinate2,q3);
J(3,4)=diff(u_coordinate2,q4);
J(3,5)=diff(u_coordinate2,q5);
J(3,6)=diff(u_coordinate2,q6);
J(3,7)=diff(u_coordinate2,x1);
J(3,8)=diff(u_coordinate2,y1);
J(3,9)=diff(u_coordinate2,z1);
J(3,10)=diff(u_coordinate2,x2);
J(3,11)=diff(u_coordinate2,y2);
J(3,12)=diff(u_coordinate2,z2);

J(4,1)=diff(v_coordinate2,q1);
J(4,2)=diff(v_coordinate2,q2);
J(4,3)=diff(v_coordinate2,q3);
J(4,4)=diff(v_coordinate2,q4);
J(4,5)=diff(v_coordinate2,q5);
J(4,6)=diff(v_coordinate2,q6);
J(4,7)=diff(v_coordinate2,x1);
J(4,8)=diff(v_coordinate2,y1);
J(4,9)=diff(v_coordinate2,z1);
J(4,10)=diff(v_coordinate2,x2);
J(4,11)=diff(v_coordinate2,y2);
J(4,12)=diff(v_coordinate2,z2);

