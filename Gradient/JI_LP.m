function [J,G_d]=JI_LPL(weight,q_cur,x_cur,Projection_cen,Im_coordinate,f_cam,d,N,M)
J=zeros(12,4,M,N);
G_d=zeros(6,M);
Gd=zeros(6,M);
Hd=zeros(6,6,M);
Hd1=zeros(3,3);
for m=1:M
    for n=1:N
        
        %u1
        J(1,1,m,n)=weight(m,n)*((x_cur(2*m-1,2)*cos(q_cur(n,1))*cos(q_cur(n,2)) + x_cur(2*m-1,3)*cos(q_cur(n,2))*sin(q_cur(n,1)))*(Im_coordinate(1,1,m,n) - Projection_cen(1,n)) - (-f_cam(n))*(x_cur(2*m-1,2)*(cos(q_cur(n,3))*sin(q_cur(n,1)) + cos(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) - x_cur(2*m-1,3)*(cos(q_cur(n,1))*cos(q_cur(n,3)) - sin(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3)))));
        J(2,1,m,n)=weight(m,n)*((-f_cam(n))*(x_cur(2*m-1,1)*sin(q_cur(n,2))*sin(q_cur(n,3)) + x_cur(2*m-1,3)*cos(q_cur(n,1))*cos(q_cur(n,2))*sin(q_cur(n,3)) - x_cur(2*m-1,2)*cos(q_cur(n,2))*sin(q_cur(n,1))*sin(q_cur(n,3))) - (Im_coordinate(1,1,m,n) - Projection_cen(1,n))*(x_cur(2*m-1,1)*cos(q_cur(n,2)) - x_cur(2*m-1,3)*cos(q_cur(n,1))*sin(q_cur(n,2)) + x_cur(2*m-1,2)*sin(q_cur(n,1))*sin(q_cur(n,2))));
        J(3,1,m,n)=weight(m,n)*(-(-f_cam(n))*(x_cur(2*m-1,2)*(cos(q_cur(n,1))*sin(q_cur(n,3)) + cos(q_cur(n,3))*sin(q_cur(n,1))*sin(q_cur(n,2))) + x_cur(2*m-1,3)*(sin(q_cur(n,1))*sin(q_cur(n,3)) - cos(q_cur(n,1))*cos(q_cur(n,3))*sin(q_cur(n,2))) + x_cur(2*m-1,1)*cos(q_cur(n,2))*cos(q_cur(n,3))));
        J(4,1,m,n)=weight(m,n)*(0);
        J(5,1,m,n)=weight(m,n)*((-f_cam(n)));
        J(6,1,m,n)=weight(m,n)*(Projection_cen(1,n) - Im_coordinate(1,1,m,n));
        J(7,1,m,n)=weight(m,n)*(- sin(q_cur(n,2))*(Im_coordinate(1,1,m,n) - Projection_cen(1,n)) - (-f_cam(n))*cos(q_cur(n,2))*sin(q_cur(n,3)));
        J(8,1,m,n)=weight(m,n)*((-f_cam(n))*(cos(q_cur(n,1))*cos(q_cur(n,3)) - sin(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) + cos(q_cur(n,2))*sin(q_cur(n,1))*(Im_coordinate(1,1,m,n) - Projection_cen(1,n)));
        J(9,1,m,n)=weight(m,n)*((-f_cam(n))*(cos(q_cur(n,3))*sin(q_cur(n,1)) + cos(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) - cos(q_cur(n,1))*cos(q_cur(n,2))*(Im_coordinate(1,1,m,n) - Projection_cen(1,n)));
        
        %v1
        J(1,2,m,n)=weight(m,n)*((-f_cam(n))*(x_cur(2*m-1,2)*(sin(q_cur(n,1))*sin(q_cur(n,3)) - cos(q_cur(n,1))*cos(q_cur(n,3))*sin(q_cur(n,2))) - x_cur(2*m-1,3)*(cos(q_cur(n,1))*sin(q_cur(n,3)) + cos(q_cur(n,3))*sin(q_cur(n,1))*sin(q_cur(n,2)))) - (x_cur(2*m-1,2)*cos(q_cur(n,1))*cos(q_cur(n,2)) + x_cur(2*m-1,3)*cos(q_cur(n,2))*sin(q_cur(n,1)))*(Projection_cen(2,n) - Im_coordinate(1,2,m,n)));
        J(2,2,m,n)=weight(m,n)*((Projection_cen(2,n) - Im_coordinate(1,2,m,n))*(x_cur(2*m-1,1)*cos(q_cur(n,2)) - x_cur(2*m-1,3)*cos(q_cur(n,1))*sin(q_cur(n,2)) + x_cur(2*m-1,2)*sin(q_cur(n,1))*sin(q_cur(n,2))) + (-f_cam(n))*(x_cur(2*m-1,1)*cos(q_cur(n,3))*sin(q_cur(n,2)) + x_cur(2*m-1,3)*cos(q_cur(n,1))*cos(q_cur(n,2))*cos(q_cur(n,3)) - x_cur(2*m-1,2)*cos(q_cur(n,2))*cos(q_cur(n,3))*sin(q_cur(n,1))));
        J(3,2,m,n)=weight(m,n)*(-(-f_cam(n))*(x_cur(2*m-1,2)*(cos(q_cur(n,1))*cos(q_cur(n,3)) - sin(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) + x_cur(2*m-1,3)*(cos(q_cur(n,3))*sin(q_cur(n,1)) + cos(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) - x_cur(2*m-1,1)*cos(q_cur(n,2))*sin(q_cur(n,3))));
        J(4,2,m,n)=weight(m,n)*(-(-f_cam(n)));
        J(5,2,m,n)=weight(m,n)*(0);
        J(6,2,m,n)=weight(m,n)*(Projection_cen(2,n) - Im_coordinate(1,2,m,n));
        J(7,2,m,n)=weight(m,n)*(sin(q_cur(n,2))*(Projection_cen(2,n) - Im_coordinate(1,2,m,n)) - (-f_cam(n))*cos(q_cur(n,2))*cos(q_cur(n,3)));
        J(8,2,m,n)=weight(m,n)*(- (-f_cam(n))*(cos(q_cur(n,1))*sin(q_cur(n,3)) + cos(q_cur(n,3))*sin(q_cur(n,1))*sin(q_cur(n,2))) - cos(q_cur(n,2))*sin(q_cur(n,1))*(Projection_cen(2,n) - Im_coordinate(1,2,m,n)));
        J(9,2,m,n)=weight(m,n)*(cos(q_cur(n,1))*cos(q_cur(n,2))*(Projection_cen(2,n) - Im_coordinate(1,2,m,n)) - (-f_cam(n))*(sin(q_cur(n,1))*sin(q_cur(n,3)) - cos(q_cur(n,1))*cos(q_cur(n,3))*sin(q_cur(n,2))));
        
        %u2
        J(1,3,m,n)=weight(m,n)*((x_cur(2*m,2)*cos(q_cur(n,1))*cos(q_cur(n,2)) + x_cur(2*m,3)*cos(q_cur(n,2))*sin(q_cur(n,1)))*(Im_coordinate(2,1,m,n) - Projection_cen(1,n)) - (-f_cam(n))*(x_cur(2*m,2)*(cos(q_cur(n,3))*sin(q_cur(n,1)) + cos(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) - x_cur(2*m,3)*(cos(q_cur(n,1))*cos(q_cur(n,3)) - sin(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3)))));
        J(2,3,m,n)=weight(m,n)*((-f_cam(n))*(x_cur(2*m,1)*sin(q_cur(n,2))*sin(q_cur(n,3)) + x_cur(2*m,3)*cos(q_cur(n,1))*cos(q_cur(n,2))*sin(q_cur(n,3)) - x_cur(2*m,2)*cos(q_cur(n,2))*sin(q_cur(n,1))*sin(q_cur(n,3))) - (Im_coordinate(2,1,m,n) - Projection_cen(1,n))*(x_cur(2*m,1)*cos(q_cur(n,2)) - x_cur(2*m,3)*cos(q_cur(n,1))*sin(q_cur(n,2)) + x_cur(2*m,2)*sin(q_cur(n,1))*sin(q_cur(n,2))));
        J(3,3,m,n)=weight(m,n)*(-(-f_cam(n))*(x_cur(2*m,2)*(cos(q_cur(n,1))*sin(q_cur(n,3)) + cos(q_cur(n,3))*sin(q_cur(n,1))*sin(q_cur(n,2))) + x_cur(2*m,3)*(sin(q_cur(n,1))*sin(q_cur(n,3)) - cos(q_cur(n,1))*cos(q_cur(n,3))*sin(q_cur(n,2))) + x_cur(2*m,1)*cos(q_cur(n,2))*cos(q_cur(n,3))));
        J(4,3,m,n)=weight(m,n)*(0);
        J(5,3,m,n)=weight(m,n)*((-f_cam(n)));
        J(6,3,m,n)=weight(m,n)*(Projection_cen(1,n) - Im_coordinate(2,1,m,n));
        J(10,3,m,n)=weight(m,n)*(- sin(q_cur(n,2))*(Im_coordinate(2,1,m,n) - Projection_cen(1,n)) - (-f_cam(n))*cos(q_cur(n,2))*sin(q_cur(n,3)));
        J(11,3,m,n)=weight(m,n)*( (-f_cam(n))*(cos(q_cur(n,1))*cos(q_cur(n,3)) - sin(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) + cos(q_cur(n,2))*sin(q_cur(n,1))*(Im_coordinate(2,1,m,n) - Projection_cen(1,n)));
        J(12,3,m,n)=weight(m,n)*((-f_cam(n))*(cos(q_cur(n,3))*sin(q_cur(n,1)) + cos(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) - cos(q_cur(n,1))*cos(q_cur(n,2))*(Im_coordinate(2,1,m,n) - Projection_cen(1,n)));
        
        %u2
        J(1,4,m,n)=weight(m,n)*((-f_cam(n))*(x_cur(2*m,2)*(sin(q_cur(n,1))*sin(q_cur(n,3)) - cos(q_cur(n,1))*cos(q_cur(n,3))*sin(q_cur(n,2))) - x_cur(2*m,3)*(cos(q_cur(n,1))*sin(q_cur(n,3)) + cos(q_cur(n,3))*sin(q_cur(n,1))*sin(q_cur(n,2)))) - (x_cur(2*m,2)*cos(q_cur(n,1))*cos(q_cur(n,2)) + x_cur(2*m,3)*cos(q_cur(n,2))*sin(q_cur(n,1)))*(Projection_cen(2,n) - Im_coordinate(2,2,m,n)));
        J(2,4,m,n)=weight(m,n)*((Projection_cen(2,n) - Im_coordinate(2,2,m,n))*(x_cur(2*m,1)*cos(q_cur(n,2)) - x_cur(2*m,3)*cos(q_cur(n,1))*sin(q_cur(n,2)) + x_cur(2*m,2)*sin(q_cur(n,1))*sin(q_cur(n,2))) + (-f_cam(n))*(x_cur(2*m,1)*cos(q_cur(n,3))*sin(q_cur(n,2)) + x_cur(2*m,3)*cos(q_cur(n,1))*cos(q_cur(n,2))*cos(q_cur(n,3)) - x_cur(2*m,2)*cos(q_cur(n,2))*cos(q_cur(n,3))*sin(q_cur(n,1))));
        J(3,4,m,n)=weight(m,n)*(-(-f_cam(n))*(x_cur(2*m,2)*(cos(q_cur(n,1))*cos(q_cur(n,3)) - sin(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) + x_cur(2*m,3)*(cos(q_cur(n,3))*sin(q_cur(n,1)) + cos(q_cur(n,1))*sin(q_cur(n,2))*sin(q_cur(n,3))) - x_cur(2*m,1)*cos(q_cur(n,2))*sin(q_cur(n,3))));
        J(4,4,m,n)=weight(m,n)*(-(-f_cam(n)));
        J(5,4,m,n)=weight(m,n)*(0);
        J(6,4,m,n)=weight(m,n)*(Projection_cen(2,n) - Im_coordinate(2,2,m,n));
        J(10,4,m,n)=weight(m,n)*(sin(q_cur(n,2))*(Projection_cen(2,n) - Im_coordinate(2,2,m,n)) - (-f_cam(n))*cos(q_cur(n,2))*cos(q_cur(n,3)));
        J(11,4,m,n)=weight(m,n)*(- (-f_cam(n))*(cos(q_cur(n,1))*sin(q_cur(n,3)) + cos(q_cur(n,3))*sin(q_cur(n,1))*sin(q_cur(n,2))) - cos(q_cur(n,2))*sin(q_cur(n,1))*(Projection_cen(2,n) - Im_coordinate(2,2,m,n)));
        J(12,4,m,n)=weight(m,n)*(cos(q_cur(n,1))*cos(q_cur(n,2))*(Projection_cen(2,n) - Im_coordinate(2,2,m,n)) - (-f_cam(n))*(sin(q_cur(n,1))*sin(q_cur(n,3)) - cos(q_cur(n,1))*cos(q_cur(n,3))*sin(q_cur(n,2))));
        
        
        
    end
    % Lagrangian constraints
    G_d(1,m)=2*x_cur(2*m-1,1) - 2*x_cur(2*m,1);
    G_d(2,m)=2*x_cur(2*m-1,2) - 2*x_cur(2*m,2);
    G_d(3,m)=2*x_cur(2*m-1,3) - 2*x_cur(2*m,3);
    G_d(4,m)=2*x_cur(2*m,1) - 2*x_cur(2*m-1,1);
    G_d(5,m)=2*x_cur(2*m,2) - 2*x_cur(2*m-1,2);
    G_d(6,m)=2*x_cur(2*m,3) - 2*x_cur(2*m-1,3);
    
    % Penalties
    % Gradient
    Gd(1,1)=+2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - d^2);
    Gd(2,m)=+2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - d^2);
    Gd(3,m)=+2*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - d^2);
    Gd(4,m)=-2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - d^2);
    Gd(5,m)=-2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - d^2);
    Gd(6,m)=-2*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3))*((x_cur(2*m-1,1) - x_cur(2*m,1))^2 + (x_cur(2*m-1,2) - x_cur(2*m,2))^2 + (x_cur(2*m-1,3) - x_cur(2*m,3))^2 - d^2);
    
    % Hessian
    Hd1(1,1)=4*(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + 4*(x_cur(2*m-1,2) - x_cur(2*m,2))^2 + 4*(x_cur(2*m-1,3) - x_cur(2*m,3))^2 + 2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))^2 - 4*d^2;
    Hd1(1,2)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2));
    Hd1(1,3)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
    
    Hd1(2,1)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2));
    Hd1(2,2)=4*(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + 4*(x_cur(2*m-1,2) - x_cur(2*m,2))^2 + 4*(x_cur(2*m-1,3) - x_cur(2*m,3))^2 + 2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))^2 - 4*d^2;
    Hd1(2,3)=2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
    
    Hd1(3,1)=2*(2*x_cur(2*m-1,1) - 2*x_cur(2*m,1))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
    Hd1(3,2)=2*(2*x_cur(2*m-1,2) - 2*x_cur(2*m,2))*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3));
    Hd1(3,3)=4*(x_cur(2*m-1,1) - x_cur(2*m,1))^2 + 4*(x_cur(2*m-1,2) - x_cur(2*m,2))^2 + 4*(x_cur(2*m-1,3) - x_cur(2*m,3))^2 + 2*(2*x_cur(2*m-1,3) - 2*x_cur(2*m,3))^2 - 4*d^2;
    
    Hd(1:6,1:6,m)=[Hd1 -Hd1;-Hd1 Hd1];
    
end

end