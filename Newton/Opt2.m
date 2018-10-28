function [Err_s,Err_f,extr_cur,wand_dis,x_cur,big_iteration]=Opt2(weight,x_cur,im_coordinate,extr_cur,P_point,focal_len,M,N_fp,N,marker_dis,itr_max,th)
%%%%%%%%%%%%%%%%%%%%%% Static Calibration Trial Run %%%%%%%%%%%%%%%%%%%%%%%
% Trial run to see which set of static parameters should be used

% -------------------------------------------------------------------------
% Very important notes :
% 1) Translation should be given in the world coordinates
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
big_iteration=0;               % Big iteration counter
delta_X=1e10;
Hd_1=zeros(3);
extr_cur1=extr_cur;
x_cur1=x_cur;


%%%%%%%%%%%%%%%%%%%%%%% Dynamic Optimization Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial value of E
Err_s=0;
for m=1:M
    for n=1:N
        temp(1)=P_point(1,n) + ((-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        temp(2)=P_point(2,n) - ((-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        temp(3)=P_point(1,n) + ((-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        temp(4)=P_point(2,n) - ((-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Err_s=Err_s+weight(m,n)*((temp(1)-im_coordinate(1,1,m,n))^2+(temp(2)-im_coordinate(1,2,m,n))^2+(temp(3)-im_coordinate(2,1,m,n))^2+(temp(4)-im_coordinate(2,2,m,n))^2);
    end
end
mu=Err_s;
P1=Err_s;
beta=1;
while (big_iteration<itr_max && max(abs(delta_X))>th)
    tic
    big_iteration=big_iteration+1
    
    %%%%%%%%%%%%%%%%% Rotation angle estimation %%%%%%%%%%%%%%%%%%
    [extr_cur]=Angle_Leven(weight,extr_cur,x_cur,N,M,im_coordinate,P_point,focal_len); % This is Levenburg-Marquardt opt. to estimate current rotation parameters
    

    %%% 3D Coordinates and translation parameters estimation %%%%%
    % Hessian matrix calculation
    Gradient=zeros(6*M+3*N,1);
    Hessian=zeros(6*M+3*N);
    for m=1:M
        for n=1:N
%             Gradient of the objective function for the m^th image of the n^th camera
            Gradient(1+(n-1)*3,1)=Gradient(1+(n-1)*3,1)+weight(m,n)*(- 2*(-focal_len(n))*((P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) - 2*(-focal_len(n))*((P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))));
            Gradient(2+(n-1)*3,1)=Gradient(2+(n-1)*3,1)+weight(m,n)*(- 2*(-focal_len(n))*((im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))) - 2*(-focal_len(n))*((im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))));
            Gradient(3+(n-1)*3,1)=Gradient(3+(n-1)*3,1)+weight(m,n)*(2*(P_point(2,n) - im_coordinate(1,2,m,n))*((P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) + 2*(P_point(2,n) - im_coordinate(2,2,m,n))*((P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) + 2*(im_coordinate(1,1,m,n) - P_point(1,n))*((im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))) + 2*(im_coordinate(2,1,m,n) - P_point(1,n))*((im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))));
            Gradient(3*N+1+(m-1)*6,1)=Gradient(3*N+1+(m-1)*6,1)+weight(m,n)*(2*((P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) + 2*((im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))*(sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3))));
            Gradient(3*N+2+(m-1)*6,1)=Gradient(3*N+2+(m-1)*6,1)+weight(m,n)*(- 2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)))*((P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) - 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)))*((im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))));
            Gradient(3*N+3+(m-1)*6,1)=Gradient(3*N+3+(m-1)*6,1)+weight(m,n)*(- 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)))*((P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) - 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)))*((im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))));
            Gradient(3*N+4+(m-1)*6,1)=Gradient(3*N+4+(m-1)*6,1)+weight(m,n)*(2*((P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) + 2*((im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))*(sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3))));
            Gradient(3*N+5+(m-1)*6,1)=Gradient(3*N+5+(m-1)*6,1)+weight(m,n)*(- 2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)))*((P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) - 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)))*((im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))));
            Gradient(3*N+6+(m-1)*6,1)=Gradient(3*N+6+(m-1)*6,1)+weight(m,n)*(- 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)))*((P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)))) - 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)))*((im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3)))));
            
%             Hessian of the objective function for the m^th image of the n^th camera
            temp1(1,1)=4*(-focal_len(n))^2;
            temp1(1,3)=- 2*(-focal_len(n))*(P_point(2,n) - im_coordinate(1,2,m,n)) - 2*(-focal_len(n))*(P_point(2,n) - im_coordinate(2,2,m,n));
            temp1(1,4)=-2*(-focal_len(n))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
            temp1(1,5)=2*(-focal_len(n))*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)));
            temp1(1,6)=2*(-focal_len(n))*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)));
            temp1(1,7)=-2*(-focal_len(n))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
            temp1(1,8)=2*(-focal_len(n))*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)));
            temp1(1,9)=2*(-focal_len(n))*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)));
            
            temp1(2,2)= 4*(-focal_len(n))^2;
            temp1(2,3)=- 2*(-focal_len(n))*(im_coordinate(1,1,m,n) - P_point(1,n)) - 2*(-focal_len(n))*(im_coordinate(2,1,m,n) - P_point(1,n));
            temp1(2,4)=-2*(-focal_len(n))*(sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            temp1(2,5)=2*(-focal_len(n))*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)));
            temp1(2,6)=2*(-focal_len(n))*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)));
            temp1(2,7)=-2*(-focal_len(n))*(sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            temp1(2,8)=2*(-focal_len(n))*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)));
            temp1(2,9)=2*(-focal_len(n))*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)));
            
            temp1(3,3)=2*(im_coordinate(1,1,m,n) - P_point(1,n))^2 + 2*(P_point(2,n) - im_coordinate(1,2,m,n))^2 + 2*(im_coordinate(2,1,m,n) - P_point(1,n))^2 + 2*(P_point(2,n) - im_coordinate(2,2,m,n))^2;
            temp1(3,4)=2*(P_point(2,n) - im_coordinate(1,2,m,n))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) + 2*(im_coordinate(1,1,m,n) - P_point(1,n))*(sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            temp1(3,5)=- 2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)))*(P_point(2,n) - im_coordinate(1,2,m,n)) - 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)))*(im_coordinate(1,1,m,n) - P_point(1,n));
            temp1(3,6)=- 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)))*(im_coordinate(1,1,m,n) - P_point(1,n)) - 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)))*(P_point(2,n) - im_coordinate(1,2,m,n));
            temp1(3,7)=2*(P_point(2,n) - im_coordinate(2,2,m,n))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) + 2*(im_coordinate(2,1,m,n) - P_point(1,n))*(sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            temp1(3,8)=- 2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)))*(P_point(2,n) - im_coordinate(2,2,m,n)) - 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)))*(im_coordinate(2,1,m,n) - P_point(1,n));
            temp1(3,9)=- 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)))*(im_coordinate(2,1,m,n) - P_point(1,n)) - 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)))*(P_point(2,n) - im_coordinate(2,2,m,n));
            
            temp1(4,4)=2*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3)))^2 + 2*(sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)))^2;
            temp1(4,5)=-2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) - 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)))*(sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            temp1(4,6)=- 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) - 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)))*(sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            
            temp1(5,5)=2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)))^2 + 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)))^2;
            temp1(5,6)=2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)))*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n))) + 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)))*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)));
            
            temp1(6,6)=2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)))^2 + 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)))^2;
            
            temp1(7,7)=2*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3)))^2 + 2*(sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)))^2;
            temp1(7,8)=- 2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) - 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)))*(sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            temp1(7,9)=- 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)))*(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3))) - 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)))*(sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) + (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            
            temp1(8,8)=2*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)))^2 + 2*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)))^2;
            temp1(8,9)=2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)))*((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n))) + 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)))*((-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)));
            
            temp1(9,9)= 2*((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)))^2 + 2*((-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)))^2;
            
            
            Hessian(1+(n-1)*3:n*3,1+(n-1)*3:n*3)=Hessian(1+(n-1)*3:n*3,1+(n-1)*3:n*3)+weight(m,n)*temp1(1:3,1:3);
            Hessian(3*N+1+(m-1)*6:3*N+m*6,3*N+1+(m-1)*6:3*N+m*6)=Hessian(3*N+1+(m-1)*6:3*N+m*6,3*N+1+(m-1)*6:3*N+m*6)+weight(m,n)*temp1(4:9,4:9);
            Hessian(1+(n-1)*3:n*3,3*N+1+(m-1)*6:3*N+m*6)=Hessian(1+(n-1)*3:n*3,3*N+1+(m-1)*6:3*N+m*6)+weight(m,n)*temp1(1:3,4:9);
     
        end
%         Gradient of the Constraint for the m^th pair of points
        Gradient(3*N+1+(m-1)*6,1)=Gradient(3*N+1+(m-1)*6,1)+2*(beta)*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2);
        Gradient(3*N+2+(m-1)*6,1)=Gradient(3*N+2+(m-1)*6,1)+2*(beta)*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2);
        Gradient(3*N+3+(m-1)*6,1)=Gradient(3*N+3+(m-1)*6,1)+2*(beta)*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2);
        Gradient(3*N+4+(m-1)*6,1)=Gradient(3*N+4+(m-1)*6,1)-2*(beta)*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2);
        Gradient(3*N+5+(m-1)*6,1)=Gradient(3*N+5+(m-1)*6,1)-2*(beta)*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2);
        Gradient(3*N+6+(m-1)*6,1)=Gradient(3*N+6+(m-1)*6,1)-2*(beta)*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2);
        
        
        % Hessian of the Constraint for the m^th pair of points
        % Hd^1_m
        Hd_1(1,1)=4*(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + 4*(x_cur(2*m-1,2) - x_cur(2*m,2))^2 + 4*(x_cur(2*m-1,3) - x_cur(2*m,3))^2 + 2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))^2 - 4*marker_dis^2;
        Hd_1(1,2)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2));
        Hd_1(1,3)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
        
        Hd_1(2,1)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2));
        Hd_1(2,2)=4*(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + 4*(x_cur(2*m-1,2) - x_cur(2*m,2))^2 + 4*(x_cur(2*m-1,3) - x_cur(2*m,3))^2 + 2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))^2 - 4*marker_dis^2;
        Hd_1(2,3)=2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
        
        Hd_1(3,1)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
        Hd_1(3,2)=2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
        Hd_1(3,3)=4*(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + 4*(x_cur(2*m-1,2) - x_cur(2*m,2))^2 + 4*(x_cur(2*m-1,3) - x_cur(2*m,3))^2 + 2*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3))^2 - 4*marker_dis^2;
        
        Hessian(3*N+1+(m-1)*6:3*N+6+(m-1)*6,3*N+1+(m-1)*6:3*N+6+(m-1)*6)=Hessian(3*N+1+(m-1)*6:3*N+6+(m-1)*6,3*N+1+(m-1)*6:3*N+6+(m-1)*6)+(beta)*[Hd_1 -Hd_1;-Hd_1 Hd_1];
    end
    Hessian=(triu(ones(6*M+3*N),+1).*Hessian)'+triu(ones(6*M+3*N)).*Hessian;
    delta_X=-(Hessian+mu*eye(6*M+3*N))\Gradient;

    % Addition of updates to the current values
    for n=1:N
        extr_cur(n,4:6)=extr_cur(n,4:6)+delta_X(3*(n-1)+1:3*n)';
    end
    for m=1:M
        x_cur(2*m-1,:)=x_cur(2*m-1,:)+delta_X(3*N+6*(m-1)+1:3*N+6*(m-1)+3)';
        x_cur(2*m,:)=x_cur(2*m,:)+delta_X(3*N+6*(m-1)+4:3*N+6*(m-1)+6)';
    end
  
    Err_f=0;
    for m=1:M
        for n=1:N
            temp(1)=P_point(1,n) + ((-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(2)=P_point(2,n) - ((-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(3)=P_point(1,n) + ((-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(4)=P_point(2,n) - ((-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            Err_f=Err_f+weight(m,n)*((temp(1)-im_coordinate(1,1,m,n))^2+(temp(2)-im_coordinate(1,2,m,n))^2+(temp(3)-im_coordinate(2,1,m,n))^2+(temp(4)-im_coordinate(2,2,m,n))^2);
        end
    end
    Err_f
    if (P1<Err_f)
        mu=2*mu;
        extr_cur=extr_cur1;
        x_cur=x_cur1;    
    elseif (P1>=Err_f)
        mu=0.5*mu;
        extr_cur1=extr_cur;
        x_cur1=x_cur;
        P1=Err_f;
        beta=1.5*beta;
    end
    toc
end


wand_dis=sqrt(sum(((x_cur(1:2:2*M,:)-x_cur(2:2:2*M,:)).^2)')'); % Marker distance for each set of 3D points.
Err_f=P1;
if (Err_f>Err_s)
    Err_f=Err_s;
end

end