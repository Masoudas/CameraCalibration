function [Err,E_im,extr_cur,wand_dis,x_cur]=Levenburg_Marquardt_LP(weight,im_coordinate,extr_cur,focal_len,P_point,M,N,marker_dis,itr_max,MSE)
%% This code optimizes the dynamic calibration function using Levenburg-Marquardt algorithm
% with Lagrange multipliers plus penalties. The object is a two marker wand wand

% -------------------------------------------------------------------------
% Very important notes :
% The projection equations that are used are:
% u=u_0-f(R(2)(x-t))/(R(3)(x-t)) ,v=v0+f(R(1)(x-t))/(R(3)(x-t))
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
itr=0;               % Big iteration counter
E_im=zeros(M,N);
Cur_im_coordinate=zeros(2,2,M,N);
for n=1:N
    extr_cur(n,4:6)=-Rotate3(extr_cur(n,1),extr_cur(n,2),extr_cur(n,3))*extr_cur(n,4:6)'; % The traslation parameter of each camera is transformed into its own coordinate system;
end
x_cur=Trg_LPL(weight,im_coordinate,extr_cur,P_point,focal_len,M,N,2); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
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
% % Distance function definition (used for penalties)
% d_fun=((x1-x2)^2+(y1-y2)^2+(z1-z2)^2-displaced_marker_dis^2)^2;
%
% % Distance function definition (used for lagrange multipliers)
% d_fun=((x1-x2)^2+(y1-y2)^2+(z1-z2)^2-displaced_marker_dis^2);
%
%
% % Cost function definition
% E=(u_coordinate1)^2+(v_coordinate1)^2+(u_coordinate2)^2+(v_coordinate2)^2;
%
% % Coordinate matrix
% C_mat=[u_coordinate1 v_coordinate1;u_coordinate2 v_coordinate2];


%%%%%%%%%%%%%%%%%%%%%%% Dynamic Optimization Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M
    for n=1:N
        Cur_im_coordinate(1,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Cur_im_coordinate(1,2,m,n)=(P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
        Cur_im_coordinate(2,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Cur_im_coordinate(2,2,m,n)=(P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
    end
end

for m=1:M
    for n=1:N
        temp(1)=P_point(1,n) - (focal_len(n)*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        temp(2)=P_point(2,n) + (focal_len(n)*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        temp(3)=P_point(1,n) - (focal_len(n)*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        temp(4)=P_point(2,n) + (focal_len(n)*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        E_im(m,n)=((im_coordinate(1,1,m,n)-temp(1))^2+(im_coordinate(1,2,m,n)-temp(2))^2)+...
            ((im_coordinate(2,1,m,n)-temp(3))^2+(im_coordinate(2,2,m,n)-temp(4))^2);
    end
end

E_im=E_im.*weight; %%% This is done since last time we decided that we do not use this image
Err=sum(sum(E_im));
Err_s=Err;
mu=Err;
% Keep the current parameters
q_cur1=extr_cur;
x_cur1=x_cur;

P1=1e18;

while ((P1-Err)>MSE && itr<itr_max)
  itr=itr+1;
    d_fun=zeros(1,M);
    
    for m=1:M
        d_fun(1,m)=(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2;
    end
   
    [J,I,Gd,Hd]=JI_LPL(weight,extr_cur,x_cur,P_point,im_coordinate,focal_len,marker_dis,N,M);
   
    %%%%%%%%%% Coefficient Matrix (H matrix) Generation %%%%%%%%%%%
    % Coefficient matrix generation
    H=zeros(6*N+6*M+M); 
    for m=1:M
        for n=1:N
            for i=1:4
                temp=J(:,i,m,n)*J(:,i,m,n)';
                H(1+6*(n-1):6*n,1+6*(n-1):6*n)=H(1+6*(n-1):6*n,1+6*(n-1):6*n)+temp(1:6,1:6);
                H(1+6*(n-1):6*n,6*N+6*(m-1)+1:6*N+6*m)=H(1+6*(n-1):6*n,6*N+6*(m-1)+1:6*N+6*m)+temp(1:6,7:12);
                H(6*N+6*(m-1)+1:6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)=H(6*N+6*(m-1)+1:6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)+temp(7:12,7:12);
                H(6*N+6*(m-1)+1:6*N+6*m,1+6*(n-1):6*n)=H(6*N+6*(m-1)+1:6*N+6*m,1+6*(n-1):6*n)+temp(7:12,1:6);
            end
        end
        H(6*N+1+6*(m-1):6*N+6*m,6*N+6*M+m)=H(6*N+1+6*(m-1):6*N+6*m,6*N+6*M+m)+1000*I(:,m);
        H(6*N+6*M+m,6*N+1+6*(m-1):6*N+6*m)=H(6*N+6*M+m,6*N+1+6*(m-1):6*N+6*m)+1000*I(:,m)';
        H(6*N+6*(m-1)+1:6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)=H(6*N+6*(m-1)+1:6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)+Hd(6,6,m);
    end
    
    % Known vector generation
    b=zeros(6*N+6*M+M,1);
    for m=1:M
        for n=1:N
            temp1=J(:,1,m,n)'*(-Cur_im_coordinate(1,1,m,n))+...
                  J(:,2,m,n)'*(-Cur_im_coordinate(1,2,m,n))+...
                  J(:,3,m,n)'*(-Cur_im_coordinate(2,1,m,n))+...
                  J(:,4,m,n)'*(-Cur_im_coordinate(2,2,m,n));
            
            b(1+6*(n-1):6*n,1)=b(1+6*(n-1):6*n,1)+weight(m,n)*temp1(1:6)';
            b(6*N+6*(m-1)+1:6*N+6*m,1)=b(6*N+6*(m-1)+1:6*N+6*m,1)+weight(m,n)*temp1(7:12)';
        end
        b(6*N+6*M+m)=-1000*d_fun(1,m);
        b(6*N+6*(m-1)+1:6*N+6*m,1)=b(6*N+6*(m-1)+1:6*N+6*m,1)+Gd(6,m);
    end
     
    H(1:6*M+6*N,1:6*M+6*N)=H(1:6*M+6*N,1:6*M+6*N)+mu*eye(6*M+6*N);
    delta_X=H\b;

    %%%%%%%%%%%%% Addition of Update to Current Values %%%%%%%%%%%%
    for n=1:N
        extr_cur(n,:)=extr_cur(n,:)+delta_X(6*(n-1)+1:6*n)';
    end
    for m=1:M
        x_cur(2*m-1,:)=x_cur(2*m-1,:)+delta_X(6*N+6*(m-1)+1:6*N+6*(m-1)+3)';
        x_cur(2*m,:)=x_cur(2*m,:)+delta_X(6*N+6*(m-1)+4:6*N+6*(m-1)+6)';
    end
    
   
    %%%%%%%%%%%%% Function Values at Each Iteration %%%%%%%%%%%%%%%
    % Updated image coordinate
    if (itr==1)
        P1=Err;
    end
    for m=1:M
        for n=1:N
            temp(1)=P_point(1,n) - (focal_len(n)*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(2)=P_point(2,n) + (focal_len(n)*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(3)=P_point(1,n) - (focal_len(n)*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(4)=P_point(2,n) + (focal_len(n)*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            E_im(m,n)=((im_coordinate(1,1,m,n)-temp(1))^2+(im_coordinate(1,2,m,n)-temp(2))^2)+...
                ((im_coordinate(2,1,m,n)-temp(3))^2+(im_coordinate(2,2,m,n)-temp(4))^2);
        end
    end
    E_im=E_im.*weight; %%% This is done since last time we decided that we do not use this image
    E2=sum(sum(E_im));
   
    %%%%%%%%% Checking for decrease in the function value %%%%%%%%%
    if (P1<E2)
        mu=2*mu;
        extr_cur=q_cur1;
        x_cur=x_cur1;
    elseif (P1>E2)
        mu=0.5*mu;
        for m=1:M
            for n=1:N
                Cur_im_coordinate(1,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
                Cur_im_coordinate(1,2,m,n)=(P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
                Cur_im_coordinate(2,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
                Cur_im_coordinate(2,2,m,n)=(P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
                E_im(m,n)=((Cur_im_coordinate(1,1,m,n))^2+(Cur_im_coordinate(1,2,m,n))^2)+...
                    ((Cur_im_coordinate(2,1,m,n))^2+(Cur_im_coordinate(2,2,m,n))^2);
            end
            
        end
        q_cur1=extr_cur;
        x_cur1=x_cur;
       P1=Err;
       Err=E2;
    end
end

wand_dis=sqrt(sum(((x_cur(1:2:2*M,:)-x_cur(2:2:2*M,:)).^2)')'); % Marker distance for each set of 3D points.
for n=1:N
    extr_cur(n,4:6)=-(Rotate3(extr_cur(n,1),extr_cur(n,2),extr_cur(n,3)))'*extr_cur(n,4:6)'; % The traslation parameter of each camera is transformed into its own coordinate system;
end



