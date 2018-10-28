function [Coordinate, im_coordinate, marker_dis]=StreoSim(N_fp,M,N,Euler_angle,translation,focal_len,P_point,Frame,Marker_displacement,Radius)

% Preprocessings ###################################
for n=1:N
    Euler_angle(n,:)=Euler_angle(n,:)*pi/180;
    Rot(1:3,1:3,n)=(Rotate3(Euler_angle(n,1),Euler_angle(n,2),Euler_angle(n,3)));
    Tr(1:3,n)=-Rot(:,:,n)*translation(n,:)';
end


% Wand #####################################
% Object Def.
Obj=[-150 0 0;150 0 0]';           % T shape

Obj=[Obj(:,1) Obj(:,2)+(Marker_displacement')];           % T shape
marker_dis=sqrt(sum((Obj(:,1)'-Obj(:,2)').^2));

Coordinate=[];
counter=0;
while (counter<=M-1)
% d_vec=(Rotate3(2*pi*rand,2*pi*rand,2*pi*rand)*[1;1;1])/sqrt(3);
% x1=[2*Radius*(rand-0.5) 2*Radius*(rand-0.5) 3000*rand];
% MObj=[x1;x1+marker_dis*d_vec']';
% if (MObj(3,2)<0)
%     MObj(3,:)=MObj(3,:)-MObj(3,2)+100*rand;
% end
if (counter<=4*M/10)
    d_vec=(Rotate3(2*pi*rand,2*pi*rand,2*pi*rand)*[1;1;1])/sqrt(3);
    x1=[500*(rand-0.5) 2*Radius*(rand-0.5) 3000*rand];
    MObj=[x1;x1+marker_dis*d_vec']';
    if (MObj(3,2)<0)
        MObj(3,:)=MObj(3,:)-MObj(3,2)+100*rand;
    end
elseif (counter>4*M/10 && counter<=4*M/5)
    d_vec=(Rotate3(2*pi*rand,2*pi*rand,2*pi*rand)*[1;1;1])/sqrt(3);
    x1=[2*Radius*(rand-0.5) 500*(rand-0.5) 3000*rand];
    MObj=[x1;x1+marker_dis*d_vec']';
    if (MObj(3,2)<0)
        MObj(3,:)=MObj(3,:)-MObj(3,2)+100*rand;
    end
    
elseif (counter>4*M/5 && counter<=M)
    d_vec=(Rotate3(2*pi*rand,2*pi*rand,2*pi*rand)*[1;1;1])/sqrt(3);
    x1=[500*(rand-0.5) 500*(rand-0.5) 3000*rand];
    MObj=[x1;x1+marker_dis*d_vec']';
    if (MObj(3,2)<0)
        MObj(3,:)=MObj(3,:)-MObj(3,2)+100*rand;
    end
end
   
 
    for n=1:N
        temp(:,:,n)=Projection(MObj,Frame(n,1:2),focal_len(n),P_point(1:2,n),Rot(1:3,1:3,n),Tr(1:3,n));
        %         ImNo = num2str(n*1000+m+1);
        %         ImName = strcat('Im',ImNo,'.bmp');
        %         imwrite(Im, ImName);
        a(n)=sum(sum(temp(:,:,n)<0));
    end

    if (sum(a==0)>=2)
       counter=counter+1;
        im_coordinate(1:N_fp,1:2,counter,1:N)=temp;
         Coordinate=[Coordinate MObj];
    end
 
end
 
Coordinate=Coordinate';

end


