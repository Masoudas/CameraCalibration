 function Err=Err_f(weight,im_coordinate,extr_cur,P_point,focal_len,M,N_fp,N)
%%%%%%%%%%%%%%%%%%%%%% Static Calibration Trial Run %%%%%%%%%%%%%%%%%%%%%%%
% Trial run to see which set of static parameters should be used

% -------------------------------------------------------------------------
% Very important notes :
% 1) Translation should be given in the world coordinates
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%% Dynamic Optimization Loop %%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial value of E
x_cur=Trg(weight,im_coordinate,extr_cur,P_point,focal_len,M,N,N_fp);
Err=0;
for m=1:M
    for n=1:N
        for j=1:N_fp
            temp(1)=P_point(1,n) + ((-focal_len(n))*(extr_cur(n,5) + x_cur(N_fp*(m-1)+j,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(N_fp*(m-1)+j,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(N_fp*(m-1)+j,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(N_fp*(m-1)+j,1)*sin(extr_cur(n,2)) + x_cur(N_fp*(m-1)+j,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(N_fp*(m-1)+j,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            temp(2)=P_point(2,n) - ((-focal_len(n))*(extr_cur(n,4) + x_cur(N_fp*(m-1)+j,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(N_fp*(m-1)+j,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(N_fp*(m-1)+j,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))))/(extr_cur(n,6) + x_cur(N_fp*(m-1)+j,1)*sin(extr_cur(n,2)) + x_cur(N_fp*(m-1)+j,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - x_cur(N_fp*(m-1)+j,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
            Err=Err+weight(m,n)*((temp(1)-im_coordinate(j,1,m,n))^2+(temp(2)-im_coordinate(j,2,m,n))^2);
        end
    end
end

end