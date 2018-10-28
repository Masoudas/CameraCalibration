function [coordinate]=Projection(Obj,FrmSiz,f,u0v0,Rot,Tr)
Mrkr_R=25;               %Marker Radius
Im=zeros(FrmSiz(1),FrmSiz(2));
[~,OP]=size(Obj);
coordinate=zeros(OP,2);

for r=1:OP
    X=[Rot Tr]*[Obj(:,r);1];
    if X(3)<0
        Imt=zeros(FrmSiz(1),FrmSiz(2));
        rc1=([(u0v0(1)-f*X(2)/X(3));(u0v0(2)+f*X(1)/X(3))]);
        rc=round([(u0v0(1)-f*X(2)/X(3));(u0v0(2)+f*X(1)/X(3))]);
        if (rc(1)>0)&&(rc(1)<FrmSiz(1))&&(rc(2)>0)&&(rc(2)<FrmSiz(2))
%             Imt(rc(1),rc(2))=1;
%             SE_R=round(max(1,Mrkr_R*f/sqrt(X'*X)));
%             SE=strel('disk',SE_R);
%             Imt=imdilate(Imt,SE);
%             Im=Im+Imt;
            coordinate(r,:)=[rc1(1) rc1(2)];
        else
            coordinate(r,:)=[-10 -10];
        end
    else
        coordinate(r,:)=-100*ones(1,2);
    end
end
end
