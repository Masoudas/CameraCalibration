function x_cur=Trg(weight,im_coordinate,extr_cur,P_point,focal_len,M,N,N_fp)
Direction_vec=zeros(N_fp,3,M,N);
temp=zeros(1,3);
x_cur=zeros(2*M,3);
for n=1:N
    Rot_mat(1:3,1:3,n)=(Rotate3(extr_cur(n,1),extr_cur(n,2),extr_cur(n,3)))'; % For each camera, the coordinate transformation matrix is generated (from the camera coordinate system to the world coordinates)
end

for i=1:N_fp
    for m=1:M
        for n=1:N
            Direction_vec(i,1:3,m,n)=(Rot_mat(:,:,n)*[im_coordinate(i,1,m,n)-P_point(1,n) im_coordinate(i,2,m,n)-P_point(2,n) focal_len(n)]'); %
            Direction_vec(i,1:3,m,n)=Direction_vec(i,1:3,m,n)/norm(Direction_vec(i,1:3,m,n));
        end
    end
end

for m=1:M % For each image
    for i=1:N_fp % For each feature point of the image
        counter=1; % Counter denotes the number of PAIR of cameras that see the point
        for j=N-1:-1:1 % From the N-1^th to the first camera
            if (weight(m,j)~=0) % If this camera sees the point
                for k=j+1:N % For each camera from j+1 to N
                    if (weight(m,k)~=0) % if this second camera sees the point as well
                        % Triangulation for the to cameras j & k
                        A(1:3,1)=Direction_vec(i,1:3,m,j)';
                        A(1:3,2)=-Direction_vec(i,1:3,m,k)';
                        b(1:3,1)=extr_cur(k,4:6)'-extr_cur(j,4:6)';
                        a=pinv(A)*b;
                        % temp(sum(0:(N-j)-1)+N-k+1,:)=a(2)*Direction_vec(i,:,m,k)+extr_cur(k,4:6);
                        temp(counter,:)=a(2)*Direction_vec(i,:,m,k)+extr_cur(k,4:6); % For this pair, save the resultant 3D location in this vecotr
                        counter=counter+1;
                        
                    end
                end
            end
        end
        for l=1:3
            x_cur((m-1)*N_fp+i,l)=mean(temp(:,l)); % Average all the 3D locations derived from triangulation.
        end
        temp=[]; % We need to empty temp to make sure that no data from the previous process affects the next estimation. This happen when we use not all images for the initial guess
    end
end
end

