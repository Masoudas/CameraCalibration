function [x_cam1,x_cam2,Dvec_world,Dvec_cam2]=StereoTrian(Im_i,Im_j,q_cur,Projection_cen,f_cam,i,j)
Dvec_world=zeros(length(Im_i),3);
Dvec_cam2=zeros(length(Im_i),3);
x_cam1=zeros(length(Im_i),3);
x_cam2=zeros(length(Im_i),3);
Rot_mat(1:3,1:3,1)=eye(3);
Rot_mat(1:3,1:3,2)=(Rotate3(q_cur(2,1),q_cur(2,2),q_cur(2,3)));
for m=1:length(Im_i)
    Dvec_world(m,1:3)=[Im_i(m,1)-Projection_cen(1,i) Im_i(m,2)-Projection_cen(2,i) f_cam(1,i)]';
    Dvec_world(m,1:3)=Dvec_world(m,1:3)/norm(Dvec_world(m,1:3));
end
for m=1:length(Im_j)
    Dvec_cam2(m,1:3)=(Rot_mat(:,:,2)'*[Im_j(m,1)-Projection_cen(1,j) Im_j(m,2)-Projection_cen(2,j) f_cam(1,j)]');
    Dvec_cam2(m,1:3)=Dvec_cam2(m,1:3)/norm(Dvec_cam2(m,1:3));
end
for m=1:length(Im_i)
    A(1:3,1)=Dvec_world(m,1:3)';
    A(1:3,2)=-Dvec_cam2(m,1:3)';
    b(1:3,1)=q_cur(2,4:6)';
    a=pinv(A)*b;
    x_cam1(m,:)=a(1)*Dvec_world(m,:,1);
    x_cam2(m,:)=(x_cam1(m,:)-q_cur(2,4:6))*Rot_mat(:,:,2)';
end
end


