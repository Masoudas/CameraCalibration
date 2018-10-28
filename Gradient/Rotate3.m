function R=Rotate3(rx,ry,rz)
Rot_x=[1 0 0; 0 cos(rx) sin(rx); 0 -sin(rx) cos(rx)];   % Rotation matrix over x axis
Rot_y=[cos(ry) 0 -sin(ry); 0 1 0; sin(ry) 0 cos(ry)];   % Rotation matrix over y axis
Rot_z=[cos(rz) sin(rz) 0; -sin(rz) cos(rz) 0; 0 0 1];   % Rotation matrix over z axis
R=Rot_z*Rot_y*Rot_x;      % Eauler rotation map
end