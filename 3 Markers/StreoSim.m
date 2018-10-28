function [im_coordinate,marker_dis]=StreoSim(N_fp,M,N,extr,Rot,focal_len,P_point,Frame,Marker_displacement,Radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates a studio by defining a moving two marker
% object with the vector 'Obj'.


% Object Def.
Obj=[-150 0 0;150 0 0]';           % T shape

Obj=[Obj(:,1) Obj(:,2)+(Marker_displacement')];           % T shape
marker_dis=sqrt(sum((Obj(:,1)'-Obj(:,2)').^2));

counter=0;
while (counter<=M-1)
    d_vec=(Rotate3(2*pi*rand,2*pi*rand,2*pi*rand)*[1;1;1])/sqrt(3);
    x1=[Radius*randn Radius*randn 3000*rand];
    MObj=[x1;x1+marker_dis*d_vec']';
    if (MObj(3,2)<0)
        MObj(3,:)=MObj(3,:)-MObj(3,2)+100*rand;
    end
   
    [temp,a]=Projection(MObj,Frame,focal_len,P_point,extr,Rot,N,N_fp);
    if (sum(a==0)>=2)
        counter=counter+1;
        im_coordinate(1:N_fp,1:2,counter,1:N)=temp;
    end
end
end


