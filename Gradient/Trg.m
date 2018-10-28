function x_cur=Trg(weight,im_coordinate,extr_cur,P_point,focal_len,M,N,N_fp)
x_cur=zeros(2*M,3);
Rot=zeros(3,3,N);
for n=1:N
    Rot(:,:,n)=Rotate3(extr_cur(n,1),extr_cur(n,2),extr_cur(n,3));
    extr_cur(n,4:6)=-(Rot(:,:,n)'*extr_cur(n,4:6)')'; % The traslation parameter of each camera is expressed in its own coordinate system;
end

for m=1:M % For each image
    for i=1:N_fp % For each feature point of the image
        t_w=0;
        temp=0;
        for j=1:N-1 % From the N-1^th to the first camera
            for k=j+1:N % For each camera from j+1 to N
                % Triangulation for the to cameras j & k
                D_vec1=(Rot(:,:,j)'*[im_coordinate(i,2,m,j)-P_point(2,j) -im_coordinate(i,1,m,j)+P_point(1,j) focal_len(j)]'); %
                D_vec2=(Rot(:,:,k)'*[im_coordinate(i,2,m,k)-P_point(2,k) -im_coordinate(i,1,m,k)+P_point(1,k) focal_len(k)]'); %
                A(1:3,1)=D_vec1'/norm(D_vec1);
                A(1:3,2)=-D_vec2'/norm(D_vec2);
                b(1:3,1)=extr_cur(k,4:6)'-extr_cur(j,4:6)';
                a=(A'*A)\(A'*b);
                temp=temp+weight(m,k)*weight(m,j)*(-a(2)*A(1:3,2)'+extr_cur(k,4:6)); % For this pair, save the resultant 3D location in this vecotr
                t_w=t_w+weight(m,j)*weight(m,k);
            end
        end
        x_cur((m-1)*N_fp+i,:)=temp/t_w; % Average all the 3D locations derived from triangulation.
    end
end
end