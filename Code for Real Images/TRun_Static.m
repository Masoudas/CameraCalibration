function [extr_1,extr_i]=TRun_Static(extr_cur,extr_s1,extr_s2,P_point,focal_len,weight,im_1i,M,N_fp,state)
% Elimination of the images not seen by the two cameras
temp1=im_1i;
im_1i=[];
counter=0;
for m=1:M    % If only one camera sees the dynamic object, we cannot use this image. So, we eliminate this image by setting its weight to zero.
    if (sum(weight(m,:))>1)
        counter=counter+1;
        weight_1i(counter,:)=weight(m,:);
        im_1i(:,:,counter,:)=temp1(:,:,m,:);
    end
end
M=counter;
if (state==2)
    Err_f11=Err_f(weight_1i,im_1i,[extr_s1(1,:,1); extr_s2(1,:,1)],P_point,focal_len,M,N_fp,2);
    Err_f12=Err_f(weight_1i,im_1i,[extr_s1(1,:,1); extr_s2(1,:,2)],P_point,focal_len,M,N_fp,2);
    Err_f21=Err_f(weight_1i,im_1i,[extr_s1(1,:,2); extr_s2(1,:,1)],P_point,focal_len,M,N_fp,2);
    Err_f22=Err_f(weight_1i,im_1i,[extr_s1(1,:,2); extr_s2(1,:,2)],P_point,focal_len,M,N_fp,2);
    temp=[Err_f11 Err_f12;Err_f21 Err_f22];
    [~,column]=min(min(temp));
    [~,row]=min(temp(:,column));
    extr_1=extr_s1(1,:,row);
    extr_i=extr_s2(1,:,column);
   
 elseif (state==1)
    Err_f11=Err_f(weight_1i,im_1i,[extr_cur(1,:); extr_s2(1,:,1)],P_point,focal_len,M,N_fp,2);
    Err_f12=Err_f(weight_1i,im_1i,[extr_cur(1,:); extr_s2(1,:,2)],P_point,focal_len,M,N_fp,2);
    if Err_f12>Err_f11
       extr_1=extr_cur(1,:);
       extr_i=extr_s2(1,:,1);
    else
       extr_1=extr_cur(1,:);
       extr_i=extr_s2(1,:,2);
    end
end


end