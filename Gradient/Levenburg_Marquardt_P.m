function [Err_s,Err_f,extr_cur,wand_dis,x_cur,itr]=Levenburg_Marquardt_P(weight,im_coordinate,extr_cur,focal_len,P_point,x_cur,M,N,marker_dis,itr_max,th)
%% This code optimizes the dynamic calibration function using Levenburg-Marquardt algorithm
% with Lagrange multipliers. The object is a two marker wand wand

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
Cur_im_coordinate=zeros(2,2,M,N);
delta_X=ones(6*M+6*N,1);


%%%%%%%%%%%%%%%%%%%%%%% Optimization Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:M
    for n=1:N
        Cur_im_coordinate(1,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Cur_im_coordinate(1,2,m,n)=(P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
        Cur_im_coordinate(2,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Cur_im_coordinate(2,2,m,n)=(P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
    end
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


Err_s=Err_f;
mu=Err_f;

% Keep the current parameters
extr_cur1=extr_cur;
x_cur1=x_cur;
P1=Err_f;


while (max(abs(delta_X(6*N+1:6*N+6*M)))>th && itr<itr_max)
    itr=itr+1;
    H=zeros(6*M+6*N);
    b=zeros(6*M+6*N,1);
    for m=1:M
        for n=1:N
            %u1
            J(1,1)=((x_cur(2*m-1,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(im_coordinate(1,1,m,n) - P_point(1,n)) - (-focal_len(n))*(x_cur(2*m-1,2)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,3)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3)))));
            J(2,1)=((-focal_len(n))*(x_cur(2*m-1,1)*sin(extr_cur(n,2))*sin(extr_cur(n,3)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*sin(extr_cur(n,3)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))*sin(extr_cur(n,3))) - (im_coordinate(1,1,m,n) - P_point(1,n))*(x_cur(2*m-1,1)*cos(extr_cur(n,2)) - x_cur(2*m-1,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m-1,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))));
            J(3,1)=(-(-focal_len(n))*(x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))));
            %         J(4,1)=(0);
            J(5,1)=((-focal_len(n)));
            J(6,1)=(P_point(1,n) - im_coordinate(1,1,m,n));
            J(7,1)=(- sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)) - (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            J(8,1)=((-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)));
            J(9,1)=((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)));
            
            %v1
            J(1,2)=((-focal_len(n))*(x_cur(2*m-1,2)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - x_cur(2*m-1,3)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2)))) - (x_cur(2*m-1,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(P_point(2,n) - im_coordinate(1,2,m,n)));
            J(2,2)=((P_point(2,n) - im_coordinate(1,2,m,n))*(x_cur(2*m-1,1)*cos(extr_cur(n,2)) - x_cur(2*m-1,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m-1,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + (-focal_len(n))*(x_cur(2*m-1,1)*cos(extr_cur(n,3))*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*cos(extr_cur(n,3)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*cos(extr_cur(n,3))*sin(extr_cur(n,1))));
            J(3,2)=(-(-focal_len(n))*(x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))));
            J(4,2)=(-(-focal_len(n)));
            %         J(5,2)=weight(m,n)*(0);
            J(6,2)=(P_point(2,n) - im_coordinate(1,2,m,n));
            J(7,2)=(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
            J(8,2)=(- (-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) - cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(1,2,m,n)));
            J(9,2)=(cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(1,2,m,n)) - (-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))));
            
            %u2
            J(1,3)=((x_cur(2*m,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(im_coordinate(2,1,m,n) - P_point(1,n)) - (-focal_len(n))*(x_cur(2*m,2)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,3)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3)))));
            J(2,3)=((-focal_len(n))*(x_cur(2*m,1)*sin(extr_cur(n,2))*sin(extr_cur(n,3)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*sin(extr_cur(n,3)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))*sin(extr_cur(n,3))) - (im_coordinate(2,1,m,n) - P_point(1,n))*(x_cur(2*m,1)*cos(extr_cur(n,2)) - x_cur(2*m,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))));
            J(3,3)=(-(-focal_len(n))*(x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))));
            %         J(4,3)=(0);
            J(5,3)=((-focal_len(n)));
            J(6,3)=(P_point(1,n) - im_coordinate(2,1,m,n));
            J(10,3)=(- sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)) - (-focal_len(n))*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
            J(11,3)=( (-focal_len(n))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)));
            J(12,3)=((-focal_len(n))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)));
            
            %u2
            J(1,4)=((-focal_len(n))*(x_cur(2*m,2)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - x_cur(2*m,3)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2)))) - (x_cur(2*m,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(P_point(2,n) - im_coordinate(2,2,m,n)));
            J(2,4)=((P_point(2,n) - im_coordinate(2,2,m,n))*(x_cur(2*m,1)*cos(extr_cur(n,2)) - x_cur(2*m,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + (-focal_len(n))*(x_cur(2*m,1)*cos(extr_cur(n,3))*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*cos(extr_cur(n,3)) - x_cur(2*m,2)*cos(extr_cur(n,2))*cos(extr_cur(n,3))*sin(extr_cur(n,1))));
            J(3,4)=(-(-focal_len(n))*(x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))));
            J(4,4)=(-(-focal_len(n)));
            %         J(5,4)=(0);
            J(6,4)=(P_point(2,n) - im_coordinate(2,2,m,n));
            J(10,4)=(sin(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
            J(11,4)=(-(-focal_len(n))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) - cos(extr_cur(n,2))*sin(extr_cur(n,1))*(P_point(2,n) - im_coordinate(2,2,m,n)));
            J(12,4)=(cos(extr_cur(n,1))*cos(extr_cur(n,2))*(P_point(2,n) - im_coordinate(2,2,m,n)) - (-focal_len(n))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))));
            
            % Hessian matrix generation
            temp1=weight(m,n)*(J(:,1)*J(:,1)'+J(:,2)*J(:,2)'+J(:,3)*J(:,3)'+J(:,4)*J(:,4)');
            
            H(1+6*(n-1):6*n,1+6*(n-1):6*n)=H(1+6*(n-1):6*n,1+6*(n-1):6*n)+temp1(1:6,1:6);
            H(1+6*(n-1):6*n,6*N+6*(m-1)+1:6*N+6*m)=H(1+6*(n-1):6*n,6*N+6*(m-1)+1:6*N+6*m)+temp1(1:6,7:12);
            H(6*N+6*(m-1)+1:6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)=H(6*N+6*(m-1)+1:6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)+temp1(7:12,7:12);
            H(6*N+6*(m-1)+1:6*N+6*m,1+6*(n-1):6*n)=H(6*N+6*(m-1)+1:6*N+6*m,1+6*(n-1):6*n)+temp1(7:12,1:6);
            
            % Known vector generation
            b(1+6*(n-1):6*n,1)=b(1+6*(n-1):6*n,1)+weight(m,n)*(J(1:6,1)*(-Cur_im_coordinate(1,1,m,n))+...
                J(1:6,2)*(-Cur_im_coordinate(1,2,m,n))+...
                J(1:6,3)*(-Cur_im_coordinate(2,1,m,n))+...
                J(1:6,4)*(-Cur_im_coordinate(2,2,m,n)));
            b(6*N+6*(m-1)+1:6*N+6*m,1)=b(6*N+6*(m-1)+1:6*N+6*m,1)+weight(m,n)*(J(7:12,1)*(-Cur_im_coordinate(1,1,m,n))+...
                J(7:12,2)*(-Cur_im_coordinate(1,2,m,n))+...
                J(7:12,3)*(-Cur_im_coordinate(2,1,m,n))+...
                J(7:12,4)*(-Cur_im_coordinate(2,2,m,n)));
        end
%         I(1,1)=2*x_cur(2*m-1,1) - 2*x_cur(2*m,1);
%         I(2,1)=2*x_cur(2*m-1,2) - 2*x_cur(2*m,2);
%         I(3,1)=2*x_cur(2*m-1,3) - 2*x_cur(2*m,3);
%         I(4,1)=2*x_cur(2*m,1) - 2*x_cur(2*m-1,1);
%         I(5,1)=2*x_cur(2*m,2) - 2*x_cur(2*m-1,2);
%         I(6,1)=2*x_cur(2*m,3) - 2*x_cur(2*m-1,3);
%         b(6*N+6*(m-1)+1:6*N+6*m,1)=b(6*N+6*(m-1)+1:6*N+6*m,1)+I*(-((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - marker_dis^2));
%         H(6*N+1+6*(m-1):6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)=H(6*N+1+6*(m-1):6*N+6*m,6*N+6*(m-1)+1:6*N+6*m)+I*I';
    end
    delta_X=(H+mu*eye(6*M+6*N))\b;
    
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

    
    %%%%%%%%% Checking for descent in the function value %%%%%%%%%
    if (P1<Err_f)
        mu=2*mu;
        extr_cur=extr_cur1;
        x_cur=x_cur1;
    elseif (P1>Err_f)
        mu=0.5*mu;
        for m=1:M
            for n=1:N
                Cur_im_coordinate(1,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(1,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
                Cur_im_coordinate(1,2,m,n)=(P_point(2,n) - im_coordinate(1,2,m,n))*(extr_cur(n,6) + x_cur(2*m-1,1)*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
                Cur_im_coordinate(2,1,m,n)=(-focal_len(n))*(extr_cur(n,5) + x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(2,1,m,n) - P_point(1,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
                Cur_im_coordinate(2,2,m,n)=(P_point(2,n) - im_coordinate(2,2,m,n))*(extr_cur(n,6) + x_cur(2*m,1)*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*(extr_cur(n,4) + x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
            end
        end
        extr_cur1=extr_cur;
        x_cur1=x_cur;
        P1=Err_f;
    end
end

if (Err_f>Err_s)
    Err_f=Err_s;
end
wand_dis=sqrt(sum(((x_cur(1:2:2*M,:)-x_cur(2:2:2*M,:)).^2)')'); % Marker distance for each set of 3D points.

end

















