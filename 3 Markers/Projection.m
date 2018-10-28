function [coordinate,a]=Projection(Obj,FrmSiz,f,u0v0,extr,Rot,N,N_fp)
coordinate=zeros(N_fp,2,N);
for n=1:N
    for r=1:N_fp
        X=[Rot(:,:,n) extr(n,4:6)']*[Obj(:,r);1];
        if X(3)<0
            rc=([(u0v0(1)-f*X(2)/X(3));(u0v0(2)+f*X(1)/X(3))]);
            if (rc(1)>0)&&(rc(1)<FrmSiz(1))&&(rc(2)>0)&&(rc(2)<FrmSiz(2))
                coordinate(r,:,n)=[rc(1) rc(2)];
            else
                coordinate(r,:,n)=[-10 -10];
            end
        else
            coordinate(r,:,n)=-100*ones(1,2);
        end
    end
    a(n)=sum(sum(coordinate(:,:,n)<0));
end
end
