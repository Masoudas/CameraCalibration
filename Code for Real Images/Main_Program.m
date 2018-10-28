
%%%%%%%%%%%%%%%%%%%% Camera Parameter Calibration %%%%%%%%%%%%%%%%%%
% This code is the main program for calibrating a multi-camera
% system using a TWO marker wand. It is assumed that all internal parameters
% are calibrated by a suitable method before this process.

% -------------------------------------------------------------------------
% Very important notes :
% 1) All length parameters are in mm.
% 2) The projection equations that are used are:
% u=u_0+f(R(1)(x-t))/(R(3)(x-t)) ,v=v0+f(R(2)(x-t))/(R(3)(x-t))
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Authors: M. Aghamohamadian-Sharbaf, H.R. Pourreza 10/6/2014
%--------------------------------------------------------------------------
close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters
N=4;                  % Number of cameras
M=300;                % Number of images (for each camera). It is assumed that all cameras take the same number of images.
N_fp=2;               % Number of feature points per image (or number of markers on the wand)
extr_s1=zeros(N,6,2); % External parameters calculated from static calibration.
MSE=1e-1;             % Desired mean square error for the cost function
mu=0;               % Newton damping parameter, in the relation (H+mu*I)
itr_s=0;            % Number of static trial run iterations
itr_max=30;        % Maximum allowable iterations for dynamic calibration optimization
marker_dis=760;

% Derived and other parameters
N_sp=N_fp*M;    % Total number of space points for all N_im images.
extr=zeros(N,6);   % This matrix contains the rotation and translation of cameras (or extrinsic parameters).
                   % Note that the first three columns for each row denote the Euler angles for its corresponding camera number.
                   % The last three parameters are the translations for that camera.
x_cur=zeros(N_sp,3);      % Current estimation of 3D points
im_coordinate=zeros(N_fp,2,M,N);   % This vector contains the coordinates of the image.
% Note that im_coordinate (:,:,m,n) denotes [u1 v1;u2 v2] of the image taken by n^th camera from the m^th set of space coordinates
extr_cur=zeros(N,6);            % This matrix contains the current estimation of extrinsic parameters parameters.
weigth=zeros(M,N);              % Weight matrix for each image of each camera. weight(m,n) denotes the weight of the m^th image for the n^th camera.
weight1=zeros(M,N);              % A template vector for weights
delta_X=zeros(3*N+3*N_fp*M,1);       % This vector denotes the update in translation and 3D parameters.
% Note that the first 3N belong to translation and the final 6M belong to space coordinates.
big_iteration=0;                % Big iteration counter
Projection_Error=1e13; % This parameter denotes the final projection error at each step. Note that this parameter is set to a nonzero value so that optimization would start.
im_static=zeros(4,2,N);
Obj_static=[0 550 0; 0 0 0; 216 0 0; 775 0 0;];
% indice=[2 3 4 5 6];
indice=[1 2 3 6];
q1=[1400:1699];
% [500:799; 800:1099; 1100:1399; 1400:1699:1700:1999]
% q=[1:q1(1)-1 1+q1(end):2200];
q2=300:499; 
% q2=1200:1299;

%%%%%%%%%%%%%%%%%%%%%%%% Camera Internal Parameters %%%%%%%%%%%%%%%%%%%%%%%
% focal_len=[1.2604e+03 1.1845e+03 1.1430e+03 1143 1.16187e+03 1.2604e+03;1.2716e+03 1.1928e+03 1.1502e+03 1169.6 1.1696e+03 1.2716e+03];  % Focal length of each camera
% P_point=[319.4175 320.0889 320.7913 320 320.7327 319.4175;240.1310 240.7918 240.6324 240 240.7711 240.1310];       % Principal point of all cameras. Note that [u_0;v_0] is the order. 
% focal_len=[ 1232.998 1120.134 1119.785 1512.391 1096.742 1145.253; 1235.102 1133.462 1118.469 1530.015 1103.427 1155.857];
% P_point=[375.790 342.842 345.652 389.956 318.420  318.816; 177.410 290.152 205.451 286.820 258.120 258.159];
% focal_len=[1.2604e+03 1.1845e+03 1.1430e+03 a 1.16187e+03 1.2604e+03;1.2716e+03 1.1928e+03 1.1502e+03 a 1.1696e+03 1.2716e+03];  % Focal length of each camera
% P_point=[319.4175 338.0889 320.7913 a 327.7327 319.4175;240.1310 308.7918 229.6324 a 243.7711 240.1310];       % Principal point of all cameras. Note that [u_0;v_0] is the order. 
focal_len=[1.2380e+03 1.1845e+03 1.1430e+03 1143 1.16187e+03 1.2604e+03;1.2503e+03 1.1928e+03 1.1502e+03 1169.6 1.1696e+03 1.2716e+03];  % Focal length of each camera
P_point=[321.2822 320.0889 320.7913 320 320.7327 319.4175;237.8460 240.7918 240.6324 240 240.7711 240.1310];       % Principal point of all cameras. Note that [u_0;v_0] is the order. 
 

%%%%%%%%%%%%%%%%%%%%%%% Loading Centers of Gravity %%%%%%%%%%%%%%%%%%%%%%%%
static=load('im_static2nd');
im_static1=static.im_static2nd;
im_static1=[im_static1(:,2,:) im_static1(:,1,:)];
n=1;
while (n<=N)
    im_static2=im_static1(:,:,indice(n));
    im_static(1:4,1:2,n)=[im_static2(:,1)-P_point(1,indice(n)) im_static2(:,2)-P_point(2,indice(n))]/focal_len(indice(n));
    % im_static(4,2,N) contains (X_c(1)/X_c(3),X_c(2)/X_c(3)) of all static images
      n=n+1;
end

dynamic=load('im_2nd');
im_coordinate1=dynamic.im;
im_coordinate1=[im_coordinate1(:,2,:,:) im_coordinate1(:,1,:,:)];
im_coordinate=im_coordinate1(:,:,q1,indice);


%%%%%%%%%%%%%%%%%%%%%%%%%  Dynamic Object Images  %%%%%%%%%%%%%%%%%%%%%%%%%
% Image elimination
% Image elimination
for m=1:M % This loop determines whether an image has less than N_fp feature points. If so, this image is eliminated
    for l=1:N
        a1=sum(sum(im_coordinate(:,:,m,l)>0));
        if (a1==N_fp*2)
            weight1(m,l)=1; % If the resultant image coordinate is negative, weight of this image is zero.
        end
    end
end

counter=0;
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


%%%%%%%%%%%%%%%%%%%%%%%%%% Static Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static calibration
for i=1:N
    [extr_s1(i,1:3,1),extr_s1(i,4:6,1),extr_s1(i,1:3,2),extr_s1(i,4:6,2)]=StaticCalib(Obj_static,im_static(:,:,i),focal_len(indice(i)));
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
    n_mat1(i,:)=max_im;
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
x_cur=Trg(weight,im_coordinate,extr_cur,P_point(:,indice),focal_len(:,indice),M,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.

%%%%%%%%%%%%%%%%%% Optimization with the Four Methods %%%%%%%%%%%%%%%%%%%%%
% Constrained Seperation method
% [Err_s_SC,Err_f_SC,extr_cur_SC,wand_dis_SC,x_cur_SC,itr_SC]=Opt_SC(weight,im_coordinate,extr_cur,P_point(:,indice),focal_len(:,indice),x_cur,M,N_fp,N,marker_dis,itr_max,MSE);
% 
% % Unconstrained Separation method
% [Err_s_S,Err_f_S,extr_cur_S,wand_dis_S,x_cur_S,itr_S]=Opt_S(weight,im_coordinate,extr_cur,P_point(:,indice),focal_len(:,indice),x_cur,M,N_fp,N,marker_dis,itr_max,MSE);
% 
% % Levenberg with Lagrange multipliers
% [Err_s_LC,Err_f_LC,extr_cur_LC,wand_dis_LC,x_cur_LC,itr_LC]=Levenburg_Marquardt_L(weight,im_coordinate,extr_cur,focal_len(:,indice),P_point(:,indice),x_cur,M,N,marker_dis,itr_max,MSE);
%                                                                                                                    
% % Levenberg with penalty functions
% [Err_s_P,Err_f_P,extr_cur_P,wand_dis_P,x_cur_P,itr_P]=Levenburg_Marquardt_P(weight,im_coordinate,extr_cur,focal_len(:,indice),P_point(:,indice),x_cur,M,N,marker_dis,itr_max,MSE);
% 
% % Unconstrained Levenberg
% [Err_s_L,Err_f_L,extr_cur_L,wand_dis_L,x_cur_L,itr_L]=Levenburg_Marquardt(weight,im_coordinate,extr_cur,focal_len(:,indice),P_point(:,indice),x_cur,M,N,marker_dis,itr_max,MSE);

n_mat=zeros(4);
for i=1:N
    for j=i+1:N
        n_mat(i,j)=sum(weight(:,j).*weight(:,i));
    end
end
n_mat
aa=input('aa=');
bb=input('bb=');
cc=input('cc=');

extr_cur_t(aa(1),1:6)=0;
extr_cur_t(aa(end),:)=StereoCalib(im_coordinate,weight,focal_len,P_point,marker_dis,aa(1),aa(2));
for n=length(bb)/2:-1:1
    if (bb(2*n-1)==bb(2*n))
        t2(n,:)=zeros(1,6);
    else
    t2(n,:)=StereoCalib(im_coordinate,weight,focal_len,P_point,marker_dis,bb(2*n-1),bb(2*n));
    end
end
Rot_temp=Rotate3(t2(length(bb)/2-1,1),t2(length(bb)/2-1,2),t2(length(bb)/2-1,3));
for o=length(bb)/2:-1:2
    Rot_temp=Rotate3(t2(o,1),t2(o,2),t2(o,3))*Rot_temp;
    tr_temp=t2(o,4:6)*Rotate3(t2(o-1,1),t2(o-1,2),t2(o-1,3))+t2(o-1,4:6);
end
x1=Rot_temp(:,:)*[1;0;0];
y1=Rot_temp(:,:)*[0;1;0];
theta(3)=phase(complex(x1(1),x1(2)));
Rot_z=[cos(-theta(3)) -sin(-theta(3)) 0; sin(-theta(3)) cos(-theta(3)) 0; 0 0 1];   % Rotating vectors around z-axis
x1=Rot_z*x1;
y1=Rot_z*y1;
theta(2)=phase(complex(x1(1),x1(3)));
Rot_y=[cos(-theta(2)) 0 sin(-theta(2)); 0 1 0; -sin(-theta(2)) 0 cos(-theta(2))];   % Rotating vectors around y-axis
y1=Rot_y'*y1;
theta(1)=phase(complex(y1(2),y1(3)));
theta(1)=-theta(1);
theta(3)=-theta(3);
extr_cur_t(bb(end),:)=[theta tr_temp];

for n=length(cc)/2:-1:1
    if (cc(2*n-1)==cc(2*n))
        t3(n,:)=zeros(1,6);
    else
    t3(n,:)=StereoCalib(im_coordinate,weight,focal_len,P_point,marker_dis,cc(2*n-1),cc(2*n));
    Rot_temp=Rotate3(t3(n,1),t3(n,2),t3(n,3));
    end
end
Rot_temp=Rotate3(t3(length(cc)/2-1,1),t3(length(cc)/2-1,2),t3(length(cc)/2-1,3));
for o=length(cc)/2:-1:2
    Rot_temp=Rotate3(t3(o,1),t3(o,2),t3(o,3))*Rot_temp;
    tr_temp=t3(o,4:6)*Rotate3(t3(o-1,1),t3(o-1,2),t3(o-1,3))+t3(o-1,4:6);
end
x1=Rot_temp(:,:)*[1;0;0];
y1=Rot_temp(:,:)*[0;1;0];
theta(3)=phase(complex(x1(1),x1(2)));
Rot_z=[cos(-theta(3)) -sin(-theta(3)) 0; sin(-theta(3)) cos(-theta(3)) 0; 0 0 1];   % Rotating vectors around z-axis
x1=Rot_z*x1;
y1=Rot_z*y1;
theta(2)=phase(complex(x1(1),x1(3)));
Rot_y=[cos(-theta(2)) 0 sin(-theta(2)); 0 1 0; -sin(-theta(2)) 0 cos(-theta(2))];   % Rotating vectors around y-axis
y1=Rot_y'*y1;
theta(1)=phase(complex(y1(2),y1(3)));
theta(1)=-theta(1);
theta(3)=-theta(3);
extr_cur_t(cc(end),:)=[theta tr_temp];

x_cur=Trg(weight,im_coordinate,extr_cur_t,P_point(:,indice),focal_len(:,indice),M,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
       pause                                                                                     % Details are given in the matlab file.
[Err_s_T,Err_f_T,extr_cur_T,wand_dis_T,x_cur_T,itr_T]=Levenburg_Marquardt_L(weight,im_coordinate,extr_cur_t,focal_len(:,indice),P_point(:,indice),x_cur,M,N,marker_dis,itr_max,MSE);

% Generating test images
M_t=100;
weight2=ones(M_t,N);
im_test=im_coordinate1(:,:,q2,indice);
for n=1:N % This loop determines whether an image has less than N_fp feature points. If so, this image is eliminated
    for m=1:M_t
        [a,~]=find(im_test(:,:,m,n)<0);
        if (~isempty(a) || length(im_test(:,1,m,n))<N_fp)
            a=[];
            weight2(m,n)=0; % If the resultant image coordinate is negative, weight of this image is zero.
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
M_t=counter2; % Number of space coordinates after elimination of useless images is given by counter.

% Static calibration initialization
x_i=Trg(weight_t,im_test,extr_cur,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_i=std(sqrt(sum(((x_i(1:2:2*M_t,:)-x_i(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_ai=mean(sqrt(sum(((x_i(1:2:2*M_t,:)-x_i(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.


% Unconstrained Seperation method
extr_cur_S(:,4:6)=extr_cur_S(:,4:6)/mean(wand_dis_S)*760;
x_s=Trg(weight_t,im_test,extr_cur_S,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_s=std(sqrt(sum(((x_s(1:2:2*M_t,:)-x_s(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_as=mean(sqrt(sum(((x_s(1:2:2*M_t,:)-x_s(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

% Constrained Seperation method
x_sc=Trg(weight_t,im_test,extr_cur_SC,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_sc=std(sqrt(sum(((x_sc(1:2:2*M_t,:)-x_sc(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_asc=mean(sqrt(sum(((x_sc(1:2:2*M_t,:)-x_sc(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

% Unconstrained Levenberg 
extr_cur_L(:,4:6)=extr_cur_L(:,4:6)/mean(wand_dis_L)*760;
x_l=Trg(weight_t,im_test,extr_cur_L,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_l=std(sqrt(sum(((x_l(1:2:2*M_t,:)-x_l(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_al=mean(sqrt(sum(((x_l(1:2:2*M_t,:)-x_l(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

% Levenberg with Lagrange multipliers
x_lc=Trg(weight_t,im_test,extr_cur_LC,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_lc=std(sqrt(sum(((x_lc(1:2:2*M_t,:)-x_lc(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_alc=mean(sqrt(sum(((x_lc(1:2:2*M_t,:)-x_lc(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

% Penalty method
x_p=Trg(weight_t,im_test,extr_cur_P,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_p=std(sqrt(sum(((x_p(1:2:2*M_t,:)-x_p(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_ap=mean(sqrt(sum(((x_p(1:2:2*M_t,:)-x_p(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

% Stereo calibration initialzation
x_t=Trg(weight_t,im_test,extr_cur_t,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_t=std(sqrt(sum(((x_t(1:2:2*M_t,:)-x_t(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_at=mean(sqrt(sum(((x_t(1:2:2*M_t,:)-x_t(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

% Levenberg method with stereo calibration initialzation
x_T=Trg(weight_t,im_test,extr_cur_T,P_point(:,indice),focal_len(:,indice),M_t,N,N_fp); % Initial 3D coordinate estimation using triangulation as explained in the thesis.
                                                                                            % Details are given in the matlab file.
wand_T=std(sqrt(sum(((x_T(1:2:2*M_t,:)-x_T(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.
wand_aT=mean(sqrt(sum(((x_T(1:2:2*M_t,:)-x_T(2:2:2*M_t,:)).^2)')')); % Marker distance for each set of 3D points.

