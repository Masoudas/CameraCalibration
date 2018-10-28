%%%%%%%%%%%%%%%%%%%% Camera Parameter Calibration %%%%%%%%%%%%%%%%%%
% This code is the main program for calibrating a multi-camera
% system using a TWO marker wand. It is assumed that all internal parameters
% are calibrated by a suitable method before this process.

% -------------------------------------------------------------------------
% Very important notes :
% 1) All length parameters are in mm.
% 2) The projection equations that are used are:
% u=u_0-f(R(2)(x-t))/(R(3)(x-t)) ,v=v0+f(R(1)(x-t))/(R(3)(x-t))
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------
close all
clear all
clc
% Loading the poses
pq1=load('pose_q1');% Random poses in the first quadrant
pq2=load('pose_q2');
pq3=load('pose_q3');
pq4=load('pose_q4');
im_s=load('im_static'); % Static images
pose(:,:,1)=pq1.pose_q1;
pose(:,:,2)=pq2.pose_q2;
pose(:,:,3)=pq3.pose_q3;
pose(:,:,4)=pq4.pose_q4;
im_static1(1:4,1:2,1:250,1:4)=im_s.im_static;
%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters
N=4;                  % Number of cameras
M=300;                % Number of images (for each camera). It is assumed that all cameras take the same number of images.
N_fp=3;               % Number of feature points per image (or number of markers on the wand)
% noise=6;              % Variance of the Gausssian noise added to all images
extr_s1=zeros(N,6,2); % External parameters calculated from static calibration.
MSE=0.1;             % Desired mean square error for the cost function
mu=0;               % Newton damping parameter, in the relation (H+mu*I)
itr_max=50;        % Maximum allowable iterations for dynamic calibration optimization
marker_displace=zeros(1,3);

% Derived and other parameters
N_sp=N_fp*M;    % Total number of space points for all N_im images.
extr=zeros(N,6);   % This matrix contains the rotation and translation of cameras (or extrinsic parameters).
% Note that the first three columns for each row denote the Euler angles for its corresponding camera number.
% The last three parameters are the translations for that camera.
im_coordinate=zeros(N_fp,2,M,N);   % This vector contains the coordinates of the image.
% Note that im_coordinate (:,:,m,n) denotes [u1 v1;u2 v2] of the image taken by n^th camera from the m^th set of space coordinates
extr_cur=zeros(N,6);            % This matrix contains the current estimation of extrinsic parameters parameters.
cam_n1=zeros(M,N);               % A template 
delta_X=zeros(3*N+3*N_fp*M,1);       % This vector denotes the update in translation and 3D parameters.
% Note that the first 3N belong to translation and the final 6M belong to space coordinates.
big_iteration=0;                % Big iteration counter
Projection_Error=1e13; % This parameter denotes the final projection error at each step. Note that this parameter is set to a nonzero value so that optimization would start.
im_static=zeros(4,2,N);
Obj_static=[0 550 0; 0 0 0; 200 0 0; 750 0 0;];
indice=1:N;

% Internal camera parameters
focal_len=984*ones(1,N);  % Focal length of each camera
P_point=[600*ones(1,N);450*ones(1,N)];       % Principal point of all cameras. Note that [u_0;v_0] is the order.
Frame=[1200*ones(N,1) 900*ones(N,1)];
for noise=[0.5 1:5]
    noise
    for test_n=1:200
        test_n      

%%%%%%%%%%%%%%%%% Camera Internal & External Parameters %%%%%%%%%%%%%%%%%%%
for k=0:N-1
    rp=randi(250,1);
    extr(k+1,1:6)=pose(rp,:,rem(k,4)+1);
    extr(k+1,1:3)=extr(k+1,1:3)/180*pi;
    n_s=noise*randn(4,2);
    im_static2=im_static1(:,:,rp,rem(k,4)+1)+n_s;
    im_static(1:4,1:2,k+1)=[im_static2(:,2)-P_point(2,k+1) -im_static2(:,1)+P_point(1,k+1)]/focal_len(k+1);
    % im_static(4,2,N) contains (X(1)/X(3),X(2)/X(3)) of all static images
    t_rad(k+1,1)=(norm(extr(k+1,4:6)));
    Rot(1:3,1:3,k+1)=(Rotate3(extr(k+1,1),extr(k+1,2),extr(k+1,3))); % For each camera, the coordinate transformation matrix is generated (from the camera coordinate system to the world coordinates)
    extr(k+1,4:6)=-(Rotate3(extr(k+1,1),extr(k+1,2),extr(k+1,3))*extr(k+1,4:6)')'; % The traslation parameter of each camera is expressed in its own coordinate system;
end
Radius=min(t_rad);

%%%%%%%%%%%%%%%%%%%%%%%%%  Dynamic Object Images  %%%%%%%%%%%%%%%%%%%%%%%%%
% Production of image coordinates
[im_coordinate, dis12, dis13, dis23]=StreoSim_3M(N_fp,M,N,extr,Rot,focal_len,P_point,Frame,marker_displace,Radius);

% Addition of noise to the images
im_coordinate=im_coordinate+noise*randn(N_fp,2,M,N);

im_temp=im_coordinate;
% Taking a subset of 200 images from the total images
M=200;
im_coordinate=im_coordinate(1:3,1:2,1:M,1:N);

weight31=zeros(300,N);              % A template vector for weights
% Image elimination
for m=1:300 % This loop determines whether an image has less than N_fp feature points. If so, this image is eliminated
    for l=1:N
        a1=sum(sum(im_temp(:,:,m,l)>0));
        if (a1==N_fp*2)
            weight31(m,l)=1; % If the resultant image coordinate is negative, weight of this image is zero.
        end
    end
end
weight3=weight31(1:M,:);

%%%%%%%%%%%%%%%%%%%%%%%%%% Static Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static calibration
for n=1:N
    [extr_s1(n,1:3,1),extr_s1(n,4:6,1),extr_s1(n,1:3,2),extr_s1(n,4:6,2)]=StaticCalib(Obj_static,im_static(:,:,n),focal_len(n));
    extr_s1(n,4:6,1)=-(Rotate3(extr_s1(n,1,1),extr_s1(n,2,1),extr_s1(n,3,1))*extr_s1(n,4:6,1)')'; % The traslation parameter of each camera is expressed in its own coordinate system;
    extr_s1(n,4:6,2)=-(Rotate3(extr_s1(n,1,2),extr_s1(n,2,2),extr_s1(n,3,2))*extr_s1(n,4:6,2)')'; % The traslation parameter of each camera is expressed in its own coordinate system;
end

% Determinig the best neighbors
n_mat=zeros(4);
for i=1:N
    max_im=zeros(1,N);
    for j=i+1:N
        max_im(j)=sum(weight3(:,j).*weight3(:,i));
    end
    max_im(i)=0;
    [n_im,index]=max(max_im);
    n_mat(i,index)=n_im;
end


% Trial-Run with static parameters
Cal_in=zeros(1,N);
while (sum(sum(n_mat))~=0)
    [~,column]=max(max(n_mat));
    [~,row]=max(n_mat(:,column));
    [~,a]=find(Cal_in==row);
    [~,b]=find(Cal_in==column);
    if (~isempty(a)&&isempty(b))
        im_1i(:,:,1:M,1)=im_coordinate(:,:,:,Cal_in(a));
        im_1i(:,:,1:M,2)=im_coordinate(:,:,:,column);
        weight_1i=[weight3(:,Cal_in(a)) weight3(:,column)];
        [~,extr_cur(column,:)]=TRun_Static(extr_cur(Cal_in(a),:),[],extr_s1(column,1:6,1:2),P_point,focal_len,weight_1i,im_1i,M,N_fp,1);
        Cal_in(column)=column;
    elseif (~isempty(b)&&isempty(a))
        im_1i(:,:,1:M,1)=im_coordinate(:,:,:,Cal_in(b));
        im_1i(:,:,1:M,2)=im_coordinate(:,:,:,row);
        weight_1i=[weight3(:,Cal_in(b)) weight3(:,row)];
        [~,extr_cur(row,:)]=TRun_Static(extr_cur(Cal_in(b),:),[],extr_s1(row,1:6,1:2),P_point,focal_len,weight_1i,im_1i,M,N_fp,1);
        Cal_in(row)=row;
    else
        im_1i(:,:,1:M,1)=im_coordinate(:,:,:,row);
        im_1i(:,:,1:M,2)=im_coordinate(:,:,:,column);
        weight_1i=[weight3(:,row) weight3(:,column)]; 
        [extr_cur(row,:),extr_cur(column,:)]=TRun_Static([],extr_s1(row,1:6,1:2),extr_s1(column,1:6,1:2),P_point,focal_len,weight_1i,im_1i,M,N_fp,2);
        Cal_in(row)=row;
        Cal_in(column)=column;
    end
    n_mat(row,column)=0;
end
warning off
%%%%%%%%%%%%%%%%%% Estimating the Initial 3D Coordinates %%%%%%%%%%%%%%%%%%
x_3M=Trg(weight3,im_coordinate,extr_cur,P_point,focal_len,M,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
                                                                     
%%%%%%%%%%%%%%%%%% Optimization with the Four Methods %%%%%%%%%%%%%%%%%%%%%
% Seperation method
[Err_s_3M(1,test_n),Err_f_3M(1,test_n),extr_cur_3M,itr_3M(ceil(noise+0.5),test_n)]=Opt_3M(weight3,im_coordinate,extr_cur,x_3M,P_point,focal_len,M,N_fp,N,dis12,dis13,dis23,itr_max,MSE);

%%%%%%%%%%%%%%%%%% Calibration with  the 2-marker object %%%%%%%%%%%%%%%%%%
M=300;
im_2marker=im_temp(1:2:3,1:2,1:M,1:N);

%%% Estimating the Initial 3D Coordinates 
x_2M=Trg(weight31,im_2marker,extr_cur,P_point,focal_len,M,N,2); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.

%%% Optimization with the Four Methods 
[Err_s_2M(1,test_n),Err_f_2M(1,test_n),extr_cur_2M,wand_dis_2M(:,test_n),x_cur_2M(:,:,test_n),itr_2M(ceil(noise+0.5),test_n)]=Opt_SC(weight31,x_2M,im_2marker,extr_cur,P_point,focal_len,M,2,N,dis13,itr_max,MSE);
                                                                                                                                               



% Generating test images
M_t=100;
[im_test,dis_t]=StreoSim(2,M_t,N,extr,Rot,focal_len,P_point,Frame,marker_displace,Radius);
                        
weightt1=zeros(M_t,N);
weightt=zeros(M_t,N);
for m=1:M_t % This loop determines whether an image has less than N_fp feature points. If so, this image is eliminated
    for n=1:N
        a1=sum(sum(im_test(:,:,m,n)>0));
        if (a1==2*2)
            weightt1(m,n)=1; % If the resultant image coordinate is negative, weight of this image is zero.
        end
    end
end

counter2=0;
temp2=im_test;
im_test=[];
for m=1:M_t    % If only one camera sees the dynamic object, we cannot use this image. So, we eliminate this image by setting its weight to zero.
    if (sum(weightt1(m,:))>1)
        counter2=counter2+1;
        weightt(counter2,:)=weightt1(m,:);
        im_test(:,:,counter2,:)=temp2(:,:,m,:);
    end
end

% Addition of noise to the images
im_test=im_test+noise*randn(2,2,M_t,N);
% No method
 
x_3M=Trg(weightt,im_test,extr_cur_3M,P_point,focal_len,M_t,N,2); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_3M(ceil(noise+0.5),test_n)=std(sqrt(sum(((x_3M(1:2:2*M_t,:)-x_3M(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
                                                                                            
wand_a3M(ceil(noise+0.5),test_n)=mean(sqrt(sum(((x_3M(1:2:2*M_t,:)-x_3M(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
% Seperation method
x_2M=Trg(weightt,im_test,extr_cur_2M,P_point,focal_len,M_t,N,2); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_2M(ceil(noise+0.5),test_n)=std(sqrt(sum(((x_2M(1:2:2*M_t,:)-x_2M(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_a2M(ceil(noise+0.5),test_n)=mean(sqrt(sum(((x_2M(1:2:2*M_t,:)-x_2M(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

end
end
