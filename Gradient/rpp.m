
function [pose,po2]=rpp(model,iprts,opt)
% Pose=rpp(model,points)	
%
% Robust Pose from Planar Tragets
% Estimates a Pose for a given Planar Target / image Points combination
% based on this 1st Solution a second solution is found.
% From both Solutions the "better one" (based on the error) is choosen
% as the correct one !
%
% (For image-sequenzes a more sophisticated approach should be used)
%
% Author: Gerald Schweighofer gerald.schweighofer@tugraz.at
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.

addpath objpose; % load path to LU
addpath util;

%% test 
% % if nargin == 0,
  model = [  0.0    0.0    200    750  ;
             550    0.0    0       0   ;
              0      0     0       0  ];

%   iprts =[  -0.0654   -0.1061   -0.1334   -0.2054;
%              0.0475    0.1117    0.1155    0.122;
%              1.0000     1.0000     1.0000     1.0000];
% end
% 
% 
% %model 3xn Model Points: planar       3rd row = 0 
% %ipts  3xn Image Points: 2d (homogen) 3rd row = 1
% 
% if nargin <= 2,
%   %% no opt -> use random values.
% %   opt.initR=rpyMat(2*pi*(rand(3,1)));
% opt.initR=[0.5502   -0.3937    0.7364
%            0.0174   -0.8763   -0.4815
%            0.8348    0.2778   -0.4753];
% end
opt.method='SVD';

%% get a first guess of the pose.
[Rlu_, tlu_, it1_, obj_err1_, img_err1_] = objpose(model, iprts(1:2,:) , opt);

%% get 2nd Pose 
sol=get2ndPose_Exact(iprts,model,Rlu_,tlu_,0);

%% refine with lu 

for i=1:length(sol),	
 opt.initR =  sol(i).R;
 [Rlu_, tlu_, it1_, obj_err1_, img_err1_] = objpose(model, iprts(1:2,:) , opt);
%  Rlu_
 sol(i).PoseLu.R = Rlu_;
 sol(i).PoseLu.t = tlu_;
 sol(i).obj_err = obj_err1_;
end

disp(['There are ' num2str(length(sol)) ' Solutions with Error: ' num2str(cat(2,sol.obj_err)) ]);
	
e = [cat(1, sol.obj_err ) [1:length(sol)]' ];
e = sortrows(e,1);

pose     =  sol(e(1,2)).PoseLu;
pose.err =  sol(e(1,2)).obj_err;

if nargout == 2,
 if size(e,1) > 1,
   po2     =  sol(e(2,2)).PoseLu;
   po2.err =  sol(e(2,2)).obj_err;
 else 
   po2 = pose;
 end
end

