function alpha=Optimum_direction_Midpoint(weight,x_cur,extr_cur,delta_X,N,M,d,beta,im_coordinate,focal_len,P_point)
interval_len=0.2;
epsilon=0.1*interval_len;
Point(1)=0;
Point(4)=1000;
E=zeros(1,4);
while (Point(4)-Point(1)>interval_len)
    mid_point=(Point(4)+Point(1))/2;
    Point(2)=mid_point-epsilon/2;
    Point(3)=mid_point+epsilon/2;
    for k=1:4
        for n=1:N
            for m=1:M
        Cur_im_coordinate(1,1)=(-focal_len(n))*((extr_cur(n,5)+Point(k)*(delta_X(3*(n-1)+2))) + (x_cur(2*m-1,2)+Point(k)*(delta_X(3*N+6*(m-1)+2)))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + (x_cur(2*m-1,3)+Point(k)*(delta_X(3*N+6*(m-1)+3)))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - (x_cur(2*m-1,1)+Point(k)*(delta_X(3*N+6*(m-1)+1)))*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(1,1,m,n) - P_point(1,n))*((extr_cur(n,6)+Point(k)*(delta_X(3*(n-1)+3))) + (x_cur(2*m-1,1)+Point(k)*(delta_X(3*N+6*(m-1)+1)))*sin(extr_cur(n,2)) + (x_cur(2*m-1,3)+Point(k)*(delta_X(3*N+6*(m-1)+3)))*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - (x_cur(2*m-1,2)+Point(k)*(delta_X(3*N+6*(m-1)+2)))*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Cur_im_coordinate(1,2)=(P_point(2,n) - im_coordinate(1,2,m,n))*((extr_cur(n,6)+Point(k)*(delta_X(3*(n-1)+3))) + (x_cur(2*m-1,1)+Point(k)*(delta_X(3*N+6*(m-1)+1)))*sin(extr_cur(n,2)) + (x_cur(2*m-1,3)+Point(k)*(delta_X(3*N+6*(m-1)+3)))*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - (x_cur(2*m-1,2)+Point(k)*(delta_X(3*N+6*(m-1)+2)))*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*((extr_cur(n,4)+Point(k)*(delta_X(3*(n-1)+1))) + (x_cur(2*m-1,2)+Point(k)*(delta_X(3*N+6*(m-1)+2)))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + (x_cur(2*m-1,3)+Point(k)*(delta_X(3*N+6*(m-1)+3)))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + (x_cur(2*m-1,1)+Point(k)*(delta_X(3*N+6*(m-1)+1)))*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
        Cur_im_coordinate(2,1)=(-focal_len(n))*((extr_cur(n,5)+Point(k)*(delta_X(3*(n-1)+2))) + (x_cur(2*m,2)+Point(k)*(delta_X(3*N+6*(m-1)+5)))*(cos(extr_cur(n,1))*cos(extr_cur(n,3)) - sin(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) + (x_cur(2*m,3)+Point(k)*(delta_X(3*N+6*(m-1)+6)))*(cos(extr_cur(n,3))*sin(extr_cur(n,1)) + cos(extr_cur(n,1))*sin(extr_cur(n,2))*sin(extr_cur(n,3))) - (x_cur(2*m,1)+Point(k)*(delta_X(3*N+6*(m-1)+4)))*cos(extr_cur(n,2))*sin(extr_cur(n,3))) - (im_coordinate(2,1,m,n) - P_point(1,n))*((extr_cur(n,6)+Point(k)*(delta_X(3*(n-1)+3))) + (x_cur(2*m,1)+Point(k)*(delta_X(3*N+6*(m-1)+4)))*sin(extr_cur(n,2)) + (x_cur(2*m,3)+Point(k)*(delta_X(3*N+6*(m-1)+6)))*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - (x_cur(2*m,2)+Point(k)*(delta_X(3*N+6*(m-1)+5)))*cos(extr_cur(n,2))*sin(extr_cur(n,1)));
        Cur_im_coordinate(2,2)=(P_point(2,n) - im_coordinate(2,2,m,n))*((extr_cur(n,6)+Point(k)*(delta_X(3*(n-1)+3))) + (x_cur(2*m,1)+Point(k)*(delta_X(3*N+6*(m-1)+4)))*sin(extr_cur(n,2)) + (x_cur(2*m,3)+Point(k)*(delta_X(3*N+6*(m-1)+6)))*cos(extr_cur(n,1))*cos(extr_cur(n,2)) - (x_cur(2*m,2)+Point(k)*(delta_X(3*N+6*(m-1)+5)))*cos(extr_cur(n,2))*sin(extr_cur(n,1))) - (-focal_len(n))*((extr_cur(n,4)+Point(k)*(delta_X(3*(n-1)+1))) + (x_cur(2*m,2)+Point(k)*(delta_X(3*N+6*(m-1)+5)))*(cos(extr_cur(n,1))*sin(extr_cur(n,3)) + cos(extr_cur(n,3))*sin(extr_cur(n,1))*sin(extr_cur(n,2))) + (x_cur(2*m,3)+Point(k)*(delta_X(3*N+6*(m-1)+6)))*(sin(extr_cur(n,1))*sin(extr_cur(n,3)) - cos(extr_cur(n,1))*cos(extr_cur(n,3))*sin(extr_cur(n,2))) + (x_cur(2*m,1)+Point(k)*(delta_X(3*N+6*(m-1)+4)))*cos(extr_cur(n,2))*cos(extr_cur(n,3)));
        E(k)=E(k)+weight(m,n)*((Cur_im_coordinate(1,1))^2+(Cur_im_coordinate(1,2))^2+(Cur_im_coordinate(2,1))^2+(Cur_im_coordinate(2,2))^2);
            end
        end
%         for m=1:M
%             E(k)=E(k)+beta*((((x_cur(2*m-1,1)+Point(k)*delta_X(3*N+6*(m-1)+1)) - (x_cur(2*m,1)+Point(k)*delta_X(3*N+6*(m-1)+4)))^2 + ((x_cur(2*m-1,2)+Point(k)*delta_X(3*N+6*(m-1)+2)) - (x_cur(2*m,2)+Point(k)*delta_X(3*N+6*(m-1)+5)))^2 + ((x_cur(2*m-1,3)+Point(k)*delta_X(3*N+6*(m-1)+3)) - (x_cur(2*m,3)+Point(k)*delta_X(3*N+6*(m-1)+6)))^2 - d^2)^2);
%         end
    end
    [~ , location(1)]=min(E);
    E(location(1))=1e16;
    [~ , location(2)]=min(E);
    if (location(1)<location(2))
        Point(1)=Point(location(1));
        Point(4)=Point(location(2));
    else
        Point(4)=Point(location(1));
        Point(1)=Point(location(2));
    end
    E=zeros(1,4);
end
alpha=(Point(1)+Point(4))/2;

end