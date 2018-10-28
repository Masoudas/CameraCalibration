function angle=Euler_ang(R)
x=R*[1;0;0];
y=R*[0;1;0];
theta_z=phase(complex(x(1),x(2)));
Rot_z=[cos(theta_z) sin(theta_z) 0; -sin(theta_z) cos(theta_z) 0; 0 0 1];   % Rotation matrix over z axis
x=Rot_z*x;
y=Rot_z*y;
theta_y=-phase(complex(x(1),x(3)));
Rot_y=[cos(theta_y) 0 -sin(theta_y); 0 1 0; sin(theta_y) 0 cos(theta_y)];   % Rotation matrix over y axis
y=Rot_y*y;
theta_x=phase(complex(y(2),y(3)));
angle=[-theta_x -theta_y -theta_z];
end