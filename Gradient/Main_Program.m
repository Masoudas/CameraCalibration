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
M=100;                % Number of images (for each camera). It is assumed that all cameras take the same number of images.
N_fp=2;               % Number of feature points per image (or number of markers on the wand)
% noise=6;              % Variance of the Gausssian noise added to all images
extr_s1=zeros(N,6,2); % External parameters calculated from static calibration.
MSE=1e-2;             % Desired mean square error for the cost function
mu=0;               % Newton damping parameter, in the relation (H+mu*I)
itr_max=300;        % Maximum allowable iterations for dynamic calibration optimization
marker_displace=zeros(1,3);

% Derived and other parameters
N_sp=N_fp*M;    % Total number of space points for all N_im images.
extr=zeros(N,6);   % This matrix contains the rotation and translation of cameras (or extrinsic parameters).
% Note that the first three columns for each row denote the Euler angles for its corresponding camera number.
% The last three parameters are the translations for that camera.
im_coordinate=zeros(N_fp,2,M,N);   % This vector contains the coordinates of the image.
% Note that im_coordinate (:,:,m,n) denotes [u1 v1;u2 v2] of the image taken by n^th camera from the m^th set of space coordinates
extr_cur=zeros(N,6);            % This matrix contains the current estimation of extrinsic parameters parameters.
weight1=zeros(M,N);              % A template vector for weights
cam_n1=zeros(M,N);               % A template 
delta_X=zeros(3*N+3*N_fp*M,1);       % This vector denotes the update in translation and 3D parameters.
% Note that the first 3N belong to translation and the final 6M belong to space coordinates.
big_iteration=0;                % Big iteration counter
counter=0;
Projection_Error=1e13; % This parameter denotes the final projection error at each step. Note that this parameter is set to a nonzero value so that optimization would start.
im_static=zeros(4,2,N);
Obj_static=[0 550 0; 0 0 0; 200 0 0; 750 0 0;];
indice=1:N;

% Internal camera parameters
focal_len=984*ones(1,N);  % Focal length of each camera
P_point=[600*ones(1,N);450*ones(1,N)];       % Principal point of all cameras. Note that [u_0;v_0] is the order.
Frame=[1200*ones(N,1) 900*ones(N,1)];
for noise=1
    noise
    for test_n=1
        test_n
      
%%%%%%%%%%%%%%%%% Camera Internal & External Parameters %%%%%%%%%%%%%%%%%%%
for n=0:N-1
    rp=randi(250,1);
    extr(n+1,1:6)=pose(rp,:,rem(n,4)+1);
    extr(n+1,1:3)=extr(n+1,1:3)/180*pi;
    n_s=noise*randn(4,2);
    im_static2=im_static1(:,:,rp,rem(n,4)+1)+n_s;
    im_static(1:4,1:2,n+1)=[im_static2(:,2)-P_point(2,n+1) -im_static2(:,1)+P_point(1,n+1)]/focal_len(n+1);
    % im_static(4,2,N) contains (X(1)/X(3),X(2)/X(3)) of all static images
    t_rad(n+1,1)=(norm(extr(n+1,4:6)));
    Rot(1:3,1:3,n+1)=(Rotate3(extr(n+1,1),extr(n+1,2),extr(n+1,3))); % For each camera, the coordinate transformation matrix is generated (from the camera coordinate system to the world coordinates)
    extr(n+1,4:6)=-(Rotate3(extr(n+1,1),extr(n+1,2),extr(n+1,3))*extr(n+1,4:6)')'; % The traslation parameter of each camera is expressed in its own coordinate system;
end
Radius=min(t_rad);

%%%%%%%%%%%%%%%%%%%%%%%%%  Dynamic Object Images  %%%%%%%%%%%%%%%%%%%%%%%%%
% Production of image coordinates
[im_coordinate,marker_dis]=StreoSim(N_fp,M,N,extr,Rot,focal_len,P_point,Frame,marker_displace,Radius);


% Image elimination
for m=1:M % This loop determines whether an image has less than N_fp feature points. If so, this image is eliminated
    for n=1:N
        a1=sum(sum(im_coordinate(:,:,m,n)>0));
        if (a1==N_fp*2)
            weight1(m,n)=1; % If the resultant image coordinate is negative, weight of this image is zero.
        end
    end
end

temp1=im_coordinate;
im_coordinate=[];
for m=1:M    % If only one camera sees the dynamic object, we cannot use this image. So, we eliminate this image by setting its weight to zero.
    if (sum(weight1(m,:))>1)
        counter=counter+1;
        weight(counter,:)=weight1(m,:);
        im_coordinate(:,:,counter,:)=temp1(:,:,m,:);
    end
end

M=counter; % Number of space coordinates after elimination of useless images is given by counter.

% Addition of noise to the images
im_coordinate=im_coordinate+noise*randn(2,2,M,N);


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
        max_im(j)=sum(weight(:,j).*weight(:,i));
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
        weight_1i=[weight(:,Cal_in(a)) weight(:,column)];
        [~,extr_cur(column,:)]=TRun_Static(extr_cur(Cal_in(a),:),[],extr_s1(column,1:6,1:2),P_point,focal_len,weight_1i,im_1i,M,N_fp,1);
        Cal_in(column)=column;
    elseif (~isempty(b)&&isempty(a))
        im_1i(:,:,1:M,1)=im_coordinate(:,:,:,Cal_in(b));
        im_1i(:,:,1:M,2)=im_coordinate(:,:,:,row);
        weight_1i=[weight(:,Cal_in(b)) weight(:,row)];
        [~,extr_cur(row,:)]=TRun_Static(extr_cur(Cal_in(b),:),[],extr_s1(row,1:6,1:2),P_point,focal_len,weight_1i,im_1i,M,N_fp,1);
        Cal_in(row)=row;
    else
        im_1i(:,:,1:M,1)=im_coordinate(:,:,:,row);
        im_1i(:,:,1:M,2)=im_coordinate(:,:,:,column);
        weight_1i=[weight(:,row) weight(:,column)]; 
        [extr_cur(row,:),extr_cur(column,:)]=TRun_Static([],extr_s1(row,1:6,1:2),extr_s1(column,1:6,1:2),P_point,focal_len,weight_1i,im_1i,M,N_fp,2);
        Cal_in(row)=row;
        Cal_in(column)=column;
    end
    n_mat(row,column)=0;
end
warning off

%%%%%%%%%%%%%%%%%% Estimating the Initial 3D Coordinates %%%%%%%%%%%%%%%%%%
x_cur=Trg(weight,im_coordinate,extr_cur,P_point,focal_len,M,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                     % Details are given in the matlab file.

                                                                     
%%%%%%%%%%%%%%%%%% Optimization with the Four Methods %%%%%%%%%%%%%%%%%%%%%
% Seperation method
tic
[Err_s_S(ceil(noise+0.5),test_n),Err_f_S(ceil(noise+0.5),test_n),extr_cur_S,wand_dis_S(:,test_n),x_cur_S(:,:,test_n),itr_S(ceil(noise+0.5),test_n)]=Opt1(weight,x_cur,im_coordinate,extr_cur,P_point,focal_len,M,N_fp,N,marker_dis,itr_max,MSE);
toc
% [Err_s_S(ceil(noise+0.5),test_n),Err_f_S(ceil(noise+0.5),test_n),extr_cur_S,wand_dis_S(:,test_n),x_cur_S(:,:,test_n),itr_S(ceil(noise+0.5),test_n)]=Opt2(weight,x_cur,im_coordinate,extr_cur,P_point,focal_len,M,N_fp,N,marker_dis,itr_max,MSE);
itr_max=100;
% [Err_s_L(ceil(noise+0.5),test_n),Err_f_L(ceil(noise+0.5),test_n),extr_cur_L,wand_dis_L(:,test_n),x_cur_L(:,:,test_n),itr_L(ceil(noise+0.5),test_n)]=Levenburg_Marquardt_L(weight,im_coordinate,extr_cur,focal_len,P_point,x_cur,M,N,marker_dis,itr_max,MSE);                                                                                                                
tic
[Err_s_P(ceil(noise+0.5),test_n),Err_f_P(ceil(noise+0.5),test_n),extr_cur_P,wand_dis_P(:,test_n),x_cur_P(:,:,test_n),itr_P(ceil(noise+0.5),test_n)]=Levenburg_Marquardt_P(weight,im_coordinate,extr_cur,focal_len,P_point,x_cur,M,N,marker_dis,itr_max,MSE);
toc


% Generating test images
M_t=100;
[im_test,dis_t]=StreoSim(N_fp,M_t,N,extr,Rot,focal_len,P_point,Frame,marker_displace,Radius);
weight2=zeros(M,N);
for m=1:M_t % This loop determines whether an image has less than N_fp feature points. If so, this image is eliminated
    for n=1:N
        a1=sum(sum(im_test(:,:,m,n)>0));
        if (a1==N_fp*2)
            weight2(m,n)=1; % If the resultant image coordinate is negative, weight of this image is zero.
        end
    end
end

counter2=0;
temp2=im_test;
im_test=[];
for m=1:M_t    % If only one camera sees the dynamic object, we cannot use this image. So, we eliminate this image by setting its weight to zero.
    if (sum(weight2(m,:))>1)
        counter2=counter2+1;
        weight_t(counter2,:)=weight2(m,:);
        im_test(:,:,counter2,:)=temp2(:,:,m,:);
    end
end

% Addition of noise to the images
im_test=im_test+noise*randn(2,2,M_t,N);
% No method
 
x_s=Trg(weight_t,im_test,extr_cur_S,P_point,focal_len,M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_s(ceil(noise+0.5),test_n)=std(sqrt(sum(((x_s(1:2:2*M_t,:)-x_s(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_as(ceil(noise+0.5),test_n)=mean(sqrt(sum(((x_s(1:2:2*M_t,:)-x_s(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
Err_s=Err_f(weight_t,im_test,extr_cur_S,P_point,focal_len,M_t,N_fp,N);
    
% Seperation method
% x_l=Trg(weight_t,im_test,extr_cur_L,P_point,focal_len,M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
%                                                                                             % Details are given in the matlab file.
% wand_l(ceil(noise+0.5),test_n)=std(sqrt(sum(((x_l(1:2:2*M_t,:)-x_l(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
% wand_al(ceil(noise+0.5),test_n)=mean(sqrt(sum(((x_l(1:2:2*M_t,:)-x_l(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
% Err_l=Err_f(weight_t,im_test,extr_cur_L,P_point,focal_len,M_t,N_fp,N);

% Seperation method
x_p=Trg(weight_t,im_test,extr_cur_P,P_point,focal_len,M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_p(ceil(noise+0.5),test_n)=std(sqrt(sum(((x_p(1:2:2*M_t,:)-x_p(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_ap(ceil(noise+0.5),test_n)=mean(sqrt(sum(((x_p(1:2:2*M_t,:)-x_p(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
Err_p=Err_f(weight_t,im_test,extr_cur_P,P_point,focal_len,M_t,N_fp,N);

end
end
