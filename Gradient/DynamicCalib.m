function [Err,E_im,extr_cur,wand_dis,x_cur]=DynamicCalib(weight,im_coordinate,extr_cur,focal_len,P_point,M,N,marker_dis,itr,MSE)
%%%%%%%%%%%%%%%%%%%%%%%% Dynamic Camera Calibration %%%%%%%%%%%%%%%%%%%%%%%
% This code is the main simulation program for calibrating a multi-camera
% system using a TWO marker wand. It is assumed that all internal parameters
% are calibrated by a suitable method before this process.

% -------------------------------------------------------------------------
% Very important notes :
% 1) All length parameters are in mm.
% 3) The projection equations that are used are:
% u=(R(1)x+t_x)/(R(3)x+t_z) ,v=(R(2)x+t_y)/(R(3)x+t_z)
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
big_iteration=0;               % Big iteration counter
E_im=zeros(M,N);
for n=1:N
    extr_cur(n,4:6)=-Rotate3(extr_cur(n,1),extr_cur(n,2),extr_cur(n,3))*extr_cur(n,4:6)'; % The traslation parameter of each camera is transformed into its own coordinate system;
end
x_cur=Trg(weight,im_coordinate,extr_cur,P_point,focal_len,M,N,2); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
% Details are given in the matlab file.


%%%%%%%%%% Basic Function Definition and Derivative Calculation %%%%%%%%%%%
% Projection function definition
% syms q1 q2 q3 q4 q5 q6 x1 y1 z1 x2 y2 z2 f0 a u_0 v_0 u01 v01 u02 v02;
% Rot_x=[1 0 0; 0 cos(q1) sin(q1); 0 -sin(q1) cos(q1)];   % Rotation matrix over x axis
% Rot_y=[cos(q2) 0 -sin(q2); 0 1 0; sin(q2) 0 cos(q2)];   % Rotation matrix over y axis
% Rot_z=[cos(q3) sin(q3) 0; -sin(q3) cos(q3) 0; 0 0 1];   % Rotation matrix over z axis
% rot_map=Rot_z*Rot_y*Rot_x;      % Eauler rotation map
% trans_vec=[q4; q5; q6];         % Translation vector of projection matrix
% Proj_mat1=[rot_map trans_vec]*[x1; y1; z1;1];
% Proj_mat2=[rot_map trans_vec]*[x2; y2; z2;1];
% u_coordinate1=f0*Proj_mat1(2,1)+(v_0-u01)*Proj_mat1(3,1);
% v_coordinate1=-f0*Proj_mat1(1,1)+(u_0-v01)*Proj_mat1(3,1);
% u_coordinate2=f0*Proj_mat2(2,1)+(v_0-u02)*Proj_mat2(3,1);
% v_coordinate2=-f0*Proj_mat2(1,1)+(u_0-v02)*Proj_mat2(3,1);
%
%
% % Distance function definition
% d_fun=((x1-x2)^2+(y1-y2)^2+(z1-z2)^2-displaced_marker_dis^2)^2;
%
% % Cost function definition
% E=(u_coordinate1)^2+(v_coordinate1)^2+(u_coordinate2)^2+(v_coordinate2)^2;
%
% % Coordinate matrix
% C_mat=[u_coordinate1 v_coordinate1;u_coordinate2 v_coordinate2];


%%%%%%%%%%%%%%%%%%%%%%% Dynamic Optimization Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
Err=0;
    for m=1:M
        for n=1:N
            temp(1)=P_point(2,n) + (focal_len(n)*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(2)=P_point(1,n) - (focal_len(n)*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(3)=P_point(2,n) + (focal_len(n)*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(4)=P_point(1,n) - (focal_len(n)*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            E_im(m,n)=weight(m,n)*((temp(1)-im_coordinate(1,1,m,n))^2+(temp(2)-im_coordinate(1,2,m,n))^2+(temp(3)-im_coordinate(2,1,m,n))^2+(temp(4)-im_coordinate(2,2,m,n))^2);
            Err=Err+weight(m,n)*((temp(1)-im_coordinate(1,1,m,n))^2+(temp(2)-im_coordinate(1,2,m,n))^2+(temp(3)-im_coordinate(2,1,m,n))^2+(temp(4)-im_coordinate(2,2,m,n))^2);
        end
    end
    Err
    mu=Err;


while (big_iteration<itr || Err>MSE)
    big_iteration=big_iteration+1;
    %%%%%%%%%%%%%%%%% Rotation angle estimation %%%%%%%%%%%%%%%%%%
    extr_cur=Angle_Leven(weight,extr_cur,x_cur,N,M,im_coordinate,P_point,focal_len); % This is Levenburg-Marquardt opt. to estimate current rotation parameters
    
    %%% 3D Coordinates and translation parameters estimation %%%%%
    % Hessian matrix calculation
    [Hessian,Gradient]=HG(weight,extr_cur,x_cur,P_point,im_coordinate,focal_len,marker_dis,N,M); % This function calculates the hessian and gradient matrix.
    
    % Update calculation
    delta_X=(Hessian+mu*eye(6*M+3*N))\(-Gradient); % Update vector calculation
    
    % Addition of updates to the current values
    for n=1:N
        extr_cur(n,4:6)=extr_cur(n,4:6)+delta_X(3*(n-1)+1:3*n)';
    end
    for m=1:M
        x_cur(2*m-1,:)=x_cur(2*m-1,:)+delta_X(3*N+6*(m-1)+1:3*N+6*(m-1)+3)';
        x_cur(2*m,:)=x_cur(2*m,:)+delta_X(3*N+6*(m-1)+4:3*N+6*(m-1)+6)';
    end
    
    %%%%%%%%%%% Projection error at each iteration %%%%%%%%%%%%%%%
    Err=0;
    for m=1:M
        for n=1:N
            temp(1)=P_point(2,n) + (focal_len(n)*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(2)=P_point(1,n) - (focal_len(n)*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(3)=P_point(2,n) + (focal_len(n)*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(4)=P_point(1,n) - (focal_len(n)*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            E_im(m,n)=weight(m,n)*((temp(1)-im_coordinate(1,1,m,n))^2+(temp(2)-im_coordinate(1,2,m,n))^2+(temp(3)-im_coordinate(2,1,m,n))^2+(temp(4)-im_coordinate(2,2,m,n))^2);
            Err=Err+weight(m,n)*((temp(1)-im_coordinate(1,1,m,n))^2+(temp(2)-im_coordinate(1,2,m,n))^2+(temp(3)-im_coordinate(2,1,m,n))^2+(temp(4)-im_coordinate(2,2,m,n))^2);
        end
    end
    Err
    mu=Err;
end


wand_dis=sqrt(sum(((x_cur(1:2:2*M,:)-x_cur(2:2:2*M,:)).^2)')'); % Marker distance for each set of 3D points.
for n=1:N
    extr_cur(n,4:6)=-(Rotate3(extr_cur(n,1),extr_cur(n,2),extr_cur(n,3)))'*extr_cur(n,4:6)'; % The traslation parameter of each camera is transformed into its own coordinate system;
end

end