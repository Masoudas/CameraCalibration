function [J,G_d]=JI_L(weight,extr_cur,x_cur,P_point,im_coordinate,focal_len,N,M)
J=zeros(12,4,M,N);
G_d=zeros(6,M);
for m=1:M
    for n=1:N
        
        %u1
        J(1,1,m,n)=weight(m,n)*((x_cur(2*m-1,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(im_coordinate(1,1,m,n) - P_point(1,n)) - focal_len(1,n)*(x_cur(2*m-1,2)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - x_cur(2*m-1,3)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2)))));
        J(2,1,m,n)=weight(m,n)*(- (im_coordinate(1,1,m,n) - P_point(1,n))*(x_cur(2*m-1,1)*cos(extr_cur(n,2)) - x_cur(2*m-1,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m-1,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))) - focal_len(1,n)*(x_cur(2*m-1,1)*cos(extr_cur(n,3))*sin(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*cos(extr_cur(n,3)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*cos(extr_cur(n,3))*sin(extr_cur(n,1))));
        J(3,1,m,n)=weight(m,n)*(focal_len(1,n)*(x_cur(2*m-1,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m-1,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))));
        J(4,1,m,n)=weight(m,n)*(focal_len(1,n));
        J(5,1,m,n)=weight(m,n)*(0);
        J(6,1,m,n)=weight(m,n)*(P_point(1,n) - im_coordinate(1,1,m,n));
        J(7,1,m,n)=weight(m,n)*(focal_len(1,n)*cos(extr_cur(n,2))*cos(extr_cur(n,3)) - sin(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)));
        J(8,1,m,n)=weight(m,n)*(focal_len(1,n)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,1,m,n) - P_point(1,n)));
        J(9,1,m,n)=weight(m,n)*(focal_len(1,n)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,1,m,n) - P_point(1,n)));

        %v1
        J(1,2,m,n)=weight(m,n)*((x_cur(2*m-1,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m-1,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(im_coordinate(1,2,m,n) - P_point(2,n)) - focal_len(2,n)*(x_cur(2*m-1,2)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m-1,3)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3)))));
        J(2,2,m,n)=weight(m,n)*(focal_len(2,n)*(x_cur(2*m-1,1)*sin(extr_cur(n,2))*sin(extr_cur(n,3)) + x_cur(2*m-1,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*sin(extr_cur(n,3)) - x_cur(2*m-1,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))*sin(extr_cur(n,3))) - (im_coordinate(1,2,m,n) - P_point(2,n))*(x_cur(2*m-1,1)*cos(extr_cur(n,2)) - x_cur(2*m-1,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m-1,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))));
        J(3,2,m,n)=weight(m,n)*(-focal_len(2,n)*(x_cur(2*m-1,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m-1,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m-1,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))));
        J(4,2,m,n)=weight(m,n)*(0);
        J(5,2,m,n)=weight(m,n)*(focal_len(2,n));
        J(6,2,m,n)=weight(m,n)*(P_point(2,n) - im_coordinate(1,2,m,n));
        J(7,2,m,n)=weight(m,n)*(- sin(extr_cur(n,2))*(im_coordinate(1,2,m,n) - P_point(2,n)) - focal_len(2,n)*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
        J(8,2,m,n)=weight(m,n)*(focal_len(2,n)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(1,2,m,n) - P_point(2,n)));
        J(9,2,m,n)=weight(m,n)*(focal_len(2,n)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(1,2,m,n) - P_point(2,n)));
        
        %u2
        J(1,3,m,n)=weight(m,n)*((x_cur(2*m,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(im_coordinate(2,1,m,n) - P_point(1,n)) - focal_len(1,n)*(x_cur(2*m,2)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - x_cur(2*m,3)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2)))));
        J(2,3,m,n)=weight(m,n)*(- (im_coordinate(2,1,m,n) - P_point(1,n))*(x_cur(2*m,1)*cos(extr_cur(n,2)) - x_cur(2*m,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))) - focal_len(1,n)*(x_cur(2*m,1)*cos(extr_cur(n,3))*sin(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*cos(extr_cur(n,3)) - x_cur(2*m,2)*cos(extr_cur(n,2))*cos(extr_cur(n,3))*sin(extr_cur(n,1))));
        J(3,3,m,n)=weight(m,n)*(focal_len(1,n)*(x_cur(2*m,2)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + x_cur(2*m,3)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,1)*cos(extr_cur(n,2))*sin(extr_cur(n,3))));
        J(4,3,m,n)=weight(m,n)*(focal_len(1,n));
        J(5,3,m,n)=weight(m,n)*(0);
        J(6,3,m,n)=weight(m,n)*(P_point(1,n) - im_coordinate(2,1,m,n));
        J(10,3,m,n)=weight(m,n)*(focal_len(1,n)*cos(extr_cur(n,2))*cos(extr_cur(n,3)) - sin(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)));
        J(11,3,m,n)=weight(m,n)*(focal_len(1,n)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,1,m,n) - P_point(1,n)));
        J(12,3,m,n)=weight(m,n)*(focal_len(1,n)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,1,m,n) - P_point(1,n)));
        
        %u2
        J(1,4,m,n)=weight(m,n)*((x_cur(2*m,2)*cos(extr_cur(n,1))*cos(extr_cur(n,2)) + x_cur(2*m,3)*cos(extr_cur(n,2))*sin(extr_cur(n,1)))*(im_coordinate(2,2,m,n) - P_point(2,n)) - focal_len(2,n)*(x_cur(2*m,2)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - x_cur(2*m,3)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3)))));
        J(2,4,m,n)=weight(m,n)*(focal_len(2,n)*(x_cur(2*m,1)*sin(extr_cur(n,2))*sin(extr_cur(n,3)) + x_cur(2*m,3)*cos(extr_cur(n,1))*cos(extr_cur(n,2))*sin(extr_cur(n,3)) - x_cur(2*m,2)*cos(extr_cur(n,2))*sin(extr_cur(n,1))*sin(extr_cur(n,3))) - (im_coordinate(2,2,m,n) - P_point(2,n))*(x_cur(2*m,1)*cos(extr_cur(n,2)) - x_cur(2*m,3)*cos(extr_cur(n,1))*sin(extr_cur(n,2)) + x_cur(2*m,2)*sin(extr_cur(n,1))*sin(extr_cur(n,2))));
        J(3,4,m,n)=weight(m,n)*(-focal_len(2,n)*(x_cur(2*m,2)*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + x_cur(2*m,3)*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + x_cur(2*m,1)*cos(extr_cur(n,2))*cos(extr_cur(n,3))));
        J(4,4,m,n)=weight(m,n)*(0);
        J(5,4,m,n)=weight(m,n)*(focal_len(2,n));
        J(6,4,m,n)=weight(m,n)*(P_point(2,n) - im_coordinate(2,2,m,n));
        J(10,4,m,n)=weight(m,n)*(- sin(extr_cur(n,2))*(im_coordinate(2,2,m,n) - P_point(2,n)) - focal_len(2,n)*cos(extr_cur(n,2))*sin(extr_cur(n,3)));
        J(11,4,m,n)=weight(m,n)*(focal_len(2,n)*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + cos(extr_cur(n,2))*sin(extr_cur(n,1))*(im_coordinate(2,2,m,n) - P_point(2,n)));
        J(12,4,m,n)=weight(m,n)*(focal_len(2,n)*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - cos(extr_cur(n,1))*cos(extr_cur(n,2))*(im_coordinate(2,2,m,n) - P_point(2,n)));
         
    end
    G_d(1,m)=2*x_cur(2*m-1,1) - 2*x_cur(2*m,1);
    G_d(2,m)=2*x_cur(2*m-1,2) - 2*x_cur(2*m,2);
    G_d(3,m)=2*x_cur(2*m-1,3) - 2*x_cur(2*m,3);
    G_d(4,m)=2*x_cur(2*m,1) - 2*x_cur(2*m-1,1);
    G_d(5,m)=2*x_cur(2*m,2) - 2*x_cur(2*m-1,2);
    G_d(6,m)=2*x_cur(2*m,3) - 2*x_cur(2*m-1,3);
    
    
end

end