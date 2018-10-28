function [angle_m,trans_m,angle2,trans2]=StaticCalib(Obj,im,f_cam)
%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%
% This function performs static calibration
% It is assumed that (X(1)/X(3),X(2)/X(3)) are given to this function and that the image plane is at z=-f;

%first estimation of depth
Lad_i=sqrt((im(1,:)-im(4,:))*(im(1,:)-im(4,:))');
Lad_Obj=sqrt((Obj(1,:)-Obj(4,:))*(Obj(1,:)-Obj(4,:))');
Dpt=Lad_Obj*f_cam/Lad_i;

%initialization
d12=sqrt((Obj(1,:)-Obj(2,:))*(Obj(1,:)-Obj(2,:))');
d13=sqrt((Obj(1,:)-Obj(3,:))*(Obj(1,:)-Obj(3,:))');
d14=sqrt((Obj(1,:)-Obj(4,:))*(Obj(1,:)-Obj(4,:))');
d23=sqrt((Obj(2,:)-Obj(3,:))*(Obj(2,:)-Obj(3,:))');
d24=sqrt((Obj(2,:)-Obj(4,:))*(Obj(2,:)-Obj(4,:))');
d34=sqrt((Obj(3,:)-Obj(4,:))*(Obj(3,:)-Obj(4,:))');

q1=[im(1,:) f_cam]; q1=q1/norm(q1);
q2=[im(2,:) f_cam]; q2=q2/norm(q2);
q3=[im(3,:) f_cam]; q3=q3/norm(q3);
q4=[im(4,:) f_cam]; q4=q4/norm(q4);

t12=q1*q2';
t13=q1*q3';
t14=q1*q4';
t23=q2*q3';
t24=q2*q4';
t34=q3*q4';


%%%%%%%%%%%%%% Optimization Loop and Parameter Estimation %%%%%%%%%%%%%%%%%
for i=-5:2
a=10^(i)*[Dpt;Dpt;Dpt;Dpt]; % Initial estimation of back-projection coef.
Da=10;
itr=1;
Lambda=10;
fap=1000*ones(6,1);
while (Da>1e-10) & itr<1000
J=2*[a(1)-t12*a(2) a(2)-t12*a(1) 0 0;...
     a(1)-t13*a(3) 0 a(3)-t13*a(1) 0;...
     a(1)-t14*a(4) 0 0 a(4)-t14*a(1);...
     0 a(2)-t23*a(3) a(3)-t23*a(2) 0;...
     0 a(2)-t24*a(4) 0 a(4)-t24*a(2);...
     0 0 a(3)-t34*a(4) a(4)-t34*a(3)];

fa=[a(1)^2-2*a(1)*a(2)*t12+a(2)^2-d12^2;...
    a(1)^2-2*a(1)*a(3)*t13+a(3)^2-d13^2;...
    a(1)^2-2*a(1)*a(4)*t14+a(4)^2-d14^2;...
    a(2)^2-2*a(2)*a(3)*t23+a(3)^2-d23^2;...
    a(2)^2-2*a(2)*a(4)*t24+a(4)^2-d24^2;...
    a(3)^2-2*a(3)*a(4)*t34+a(4)^2-d34^2];
Delta=-1*inv(J'*J+Lambda*eye(4,4))*J'*fa;

a=a+Delta;
itr=itr+1;
    if min(abs(fa))<min(abs(fap))
    Lambda=10*Lambda;
else
    Lambda=Lambda/10;
end
fap=fa;
end
A=q1*a(1);
B=q2*a(2);
C=q3*a(3);
D=q4*a(4);

AB=sqrt((A-B)*(A-B)');
BC=sqrt((B-C)*(B-C)');
BD=sqrt((B-D)*(B-D)');
CD=sqrt((C-D)*(C-D)');
AD=sqrt((A-D)*(A-D)');

dis1(:,i+6)=sum([(AB-550)^2 (BC-200)^2 (BD-750)^ (CD-550)^2 (AD-sqrt(550^2+750^2))^2]); % Sum of the difference between estimated and actual length

    %%% External parameter estimation for current initial guess
T=-B';
y=(A-B)'/norm((A-B)); % y axis
x=(D-B)'/norm(D-B); % x axis

theta_z=angle(complex(x(1),x(2)));
Rot_z=[cos(theta_z) sin(theta_z) 0; -sin(theta_z) cos(theta_z) 0; 0 0 1];   % Rotation matrix over z axis
x=Rot_z*x;
y=Rot_z*y;
theta_y=-angle(complex(x(1),x(3)));
Rot_y=[cos(theta_y) 0 -sin(theta_y); 0 1 0; sin(theta_y) 0 cos(theta_y)];   % Rotation matrix over y axis
y=Rot_y*y;
theta_x=angle(complex(y(2),y(3)));
Rot_x=[1 0 0; 0 cos(theta_x) sin(theta_x); 0 -sin(theta_x) cos(theta_x)];   % Rotation matrix over x axis
tet(i+6,:)=[-theta_x -theta_y -theta_z]/pi*180;
Translation1(:,i+6)=(Rot_x*Rot_y*Rot_z*(T));

end
j=1;
for i=1:8
    if (norm(Translation1(:,i))<10000) % If a camera is estimated to be more than 10m away from the Jig, the estimation is incorrect
        Translation(:,j)=Translation1(:,i);
        angle1(j,:)=tet(i,:);
        dis(:,j)=dis1(:,i);
        j=j+1;
    end
end


[~,j]=min(dis);
angle1=angle1(j,:);
Translation=Translation(:,j);
iprts=[im';ones(1,4)];
R=Rotate3(angle1(1)/180*pi,angle1(2)/180*pi,angle1(3)/180*pi);
opt=R;

% Parameters of the main guess
[pose,po2]=rpp(Obj,iprts,opt);
R_m=[-pose.R(:,1) -pose.R(:,2) pose.R(:,3)];
trans_m=R_m'*pose.t;
angle_m=Euler_ang(R_m);

% Parameters of the second guess
R_2=[-po2.R(:,1) -po2.R(:,2) po2.R(:,3)];
trans2=R_2'*po2.t;
angle2=Euler_ang(R_2);

end