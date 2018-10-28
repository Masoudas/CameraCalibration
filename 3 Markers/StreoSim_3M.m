function [im_coordinate, dis12, dis13, dis23]=StreoSim_3M(N_fp,M,N,extr,Rot,focal_len,P_point,Frame,Marker_displacement,Radius)

% Wand #####################################
% Object Def.
Obj=[-150 0 0;50 0 0;150 0 0]';           % T shape
dis12=sqrt(sum((Obj(:,1)'-Obj(:,2)').^2));
dis13=sqrt(sum((Obj(:,1)'-Obj(:,3)').^2));
dis23=sqrt(sum((Obj(:,2)'-Obj(:,3)').^2));


counter=0;
while (counter<=M-1)
        d_vec=(Rotate3(2*pi*rand,2*pi*rand,2*pi*rand)*[1;1;1])/sqrt(3);
        x1=[500*(rand-0.5) 2*Radius*(rand-0.5) 3000*rand];
        MObj=[x1;x1+dis12*d_vec';x1+dis13*d_vec']';
        if (sum(MObj(3,:)<0)>1)
            MObj(3,:)=MObj(3,:)-MObj(3,2)+100*rand;
        end
    [temp,a]=Projection(MObj,Frame,focal_len,P_point,extr,Rot,N,N_fp);
    if (sum(a==0)>=2)
        counter=counter+1;
        im_coordinate(1:N_fp,1:2,counter,1:N)=temp;
    end
end
    

end


