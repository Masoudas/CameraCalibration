function [extr]=StereoCalib(im_coordinate,weight,focal_len,P_point,Marker_distance,index1,index2)
% This code is calibrates two cameras with respect to each other
% using a TWO marker wand. It is assumed that all internal parameters
% are calibrated by a suitable method before this process.
% Very important note:
% This function works for images with two markers

% -------------------------------------------------------------------------
% Very important notes :
% u=u0+f*X(1)/X(3) , v=v0+f*X(2)/X(3)
% 
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------


% Finding Matching Points 
k=1;
for m=1:length(weight(:,index1))
    if weight(m,index1)>0 && weight(m,index2)>0
        for n=1:2
            Im_1(2*(k-1)+n,:)=im_coordinate(n,:,m,index1);
            Im_2(2*(k-1)+n,:)=im_coordinate(n,:,m,index2);
        end
         k=k+1;
    end
end

% Estimating fundamental matrix
[F,index,status]=estimateFundamentalMatrix(Im_1,Im_2,'Method', 'RANSAC', 'DistanceThreshold', 1);

if (status==0) % If F is estimated
K1=[focal_len(index1) 0 P_point(1,index1);0 focal_len(index1) P_point(2,index1); 0 0 1];
K2=[focal_len(index2) 0 P_point(1,index2);0 focal_len(index2) P_point(2,index2); 0 0 1];
Essential_mat=K2'*F*K1;
E=[0 1 0;-1 0 0;0 0 1];
[U,S,V]=svd(Essential_mat);

if (det(U*V')<0) % if singular values are negative, change their sign
    [U,S,V]=svd(-Essential_mat);
end

% Determining four sets of parameters
R(:,:,1)=(U*E*V');
R(:,:,2)=(U*E'*V');
R(:,:,3)=(U*E*V');
R(:,:,4)=(U*E'*V');
Translation_prime(:,1)=-R(:,:,1)'*U*[0;0;1];
Translation_prime(:,2)=-R(:,:,2)'*U*[0;0;1];
Translation_prime(:,3)=R(:,:,3)'*U*[0;0;1];
Translation_prime(:,4)=R(:,:,4)'*U*[0;0;1];
theta=zeros(4,3);

% Finding the rotation angle and 3D coordinates for each set
for l=1:4
    x1=R(:,:,l)*[1;0;0];
    y1=R(:,:,l)*[0;1;0];
    theta(l,3)=phase(complex(x1(1),x1(2)));
    Rot_z=[cos(-theta(l,3)) -sin(-theta(l,3)) 0; sin(-theta(l,3)) cos(-theta(l,3)) 0; 0 0 1];   % Rotating vectors around z-axis
    x1=Rot_z*x1;
    y1=Rot_z*y1;
    theta(l,2)=phase(complex(x1(1),x1(3)));
    Rot_y=[cos(-theta(l,2)) 0 sin(-theta(l,2)); 0 1 0; -sin(-theta(l,2)) 0 cos(-theta(l,2))];   % Rotating vectors around y-axis
    y1=Rot_y'*y1;
    theta(l,1)=phase(complex(y1(2),y1(3)));
    theta(l,1)=-theta(l,1);
    theta(l,3)=-theta(l,3);
    [x_cam1(:,:,l),x_cam2(:,:,l)]=StereoTrian(Im_1,Im_2,[zeros(1,6);theta(l,1) theta(l,2) theta(l,3) Translation_prime(:,l)'],P_point,focal_len,index1,index2);
end
% Determinig the correct set of parameters 
for l=1:4
    count(l)=sum(x_cam1(index,3,l)>0)+sum(x_cam2(index,3,l)>0);
end
[~,c_index]=max(count);

% Determining the Scale of the Scene 
for l=1:4
    X(:,l)=(sqrt(sum(((x_cam1(1:2:length(x_cam1(:,1,1)),:,l)-x_cam1(2:2:length(x_cam1(:,1,1)),:,l)).*(x_cam1(1:2:length(x_cam1(:,1,1)),:,l)-x_cam1(2:2:length(x_cam1(:,1,1)),:,l)))')))';
end
for l=1:4
    Translation(:,l)=Marker_distance*Translation_prime(:,l)/(mean(X(:,l)));
end

rot=theta(c_index,:);
trans=Translation(:,c_index);
extr=[rot trans'];

else
    print('Calibration is impossible')
end
end