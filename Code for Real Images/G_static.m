% Noghte 1, 0 550 0, 2 0 0 0, 3 200 0 0, ...
clear all
clc
for counter=1:6
    
    f=imread(strcat('E:\Sequences\LFrame-0',num2str(counter),'-03.jpg'));
    g=rgb2gray(f);
    [R C]=size(g);
    level = graythresh(g);
    gb = im2bw(g, 0.5);
    [L,BP] = bwlabel(gb);
    COG=[];
    Area=[];
    for i=1:BP
        L1=(L==i);
        A=sum(L1(:));
        x0=0; y0=0;
        for r=1:R
            for c=1:C
                x0=x0+c*L1(r,c);
                y0=y0+r*L1(r,c);
            end
        end
        x0=x0/A; y0=y0/A;
        COG=[COG;x0 y0];
        Area=[Area; A];
    end
    COGt=COG;
    [As,I] = sort(Area, 'descend');
    COG=[];
    for r=1:4
        COG(r,:)=COGt(I(r),:);
    end
    
    % ##############################################################
    % FindABCD
    % Find A
    ABCD=zeros(4,2);
    CL=[];
    a=COG(1,:);b=COG(2,:);c=COG(3,:);
    CL=[CL ;(c(1)-a(1))*(b(2)-a(2))+(c(2)-a(2))*(a(1)-b(1))];
    a=COG(1,:);b=COG(2,:);c=COG(4,:);
    CL=[CL ;(c(1)-a(1))*(b(2)-a(2))+(c(2)-a(2))*(a(1)-b(1))];
    a=COG(1,:);b=COG(3,:);c=COG(4,:);
    CL=[CL ;(c(1)-a(1))*(b(2)-a(2))+(c(2)-a(2))*(a(1)-b(1))];
    a=COG(2,:);b=COG(3,:);c=COG(4,:);
    CL=[CL ;(c(1)-a(1))*(b(2)-a(2))+(c(2)-a(2))*(a(1)-b(1))];
    CL=abs(CL);
    [CLMin I]=min(CL);
    switch I
        case 1.0
            ABCD(1,:)=COG(4,:); a1=1; b1=2; c1=3;
        case 2.0
            ABCD(1,:)=COG(3,:); a1=1; b1=2; c1=4;
        case 3.0
            ABCD(1,:)=COG(2,:); a1=1; b1=3; c1=4;
        otherwise
            ABCD(1,:)=COG(1,:); a1=2; b1=3; c1=4;
    end
    
%     % Find D
%     D=[(ABCD(1,:)-COG(a1,:))*(ABCD(1,:)-COG(a1,:))'; ...
%         (ABCD(1,:)-COG(b1,:))*(ABCD(1,:)-COG(b1,:))'; ...
%         (ABCD(1,:)-COG(c1,:))*(ABCD(1,:)-COG(c1,:))'];
%     [DMax I]=max(D);
%     switch I
%         case 1.0
%             ABCD(4,:)=COG(a1,:); a2=b1; b2=c1;
%         case 2.0
%             ABCD(4,:)=COG(b1,:); a2=a1; b2=c1;
%         otherwise
%             ABCD(4,:)=COG(c1,:); a2=a1; b2=b1;
%     end
%     
%     % Find BC
%     D=[(ABCD(4,:)-COG(a2,:))*(ABCD(4,:)-COG(a2,:))'; ...
%         (ABCD(4,:)-COG(b2,:))*(ABCD(4,:)-COG(b2,:))'];
%     [DMax I]=max(D);
%     switch I
%         case 1.0
%             ABCD(2,:)=COG(a2,:); ABCD(3,:)=COG(b2,:);
%         otherwise
%             ABCD(2,:)=COG(b2,:); ABCD(3,:)=COG(a2,:);
%     end
     D=zeros(2,2);
        dis1(1)=norm(COG(a1,:)-COG(b1,:));
        dis1(2)=norm(COG(a1,:)-COG(c1,:));
        dis1(3)=norm(COG(b1,:)-COG(c1,:));
        [~, i1]=max(dis1);
        if (i1==1)
            if (dis1(2)>dis1(3))
                ABCD(4,1:2)=COG(a1,:);
                ABCD(2,1:2)=COG(b1,:);
                ABCD(3,1:2)=COG(c1,:);
            elseif (dis1(2)<dis1(3))
                ABCD(4,1:2)=COG(b1,:);
                ABCD(2,1:2)=COG(a1,:);
                ABCD(3,1:2)=COG(c1,:);
            end
        elseif (i1==2)
            if (dis1(1)>dis1(3))
                ABCD(4,1:2)=COG(a1,:);
                ABCD(2,1:2)=COG(c1,:);
                ABCD(3,1:2)=COG(b1,:);
            elseif (dis1(1)<dis1(3))
                ABCD(4,1:2)=COG(c1,:);
                ABCD(2,1:2)=COG(a1,:);
                ABCD(3,1:2)=COG(b1,:);
            end
        elseif (i1==3)
            if (dis1(1)>dis1(2))
                ABCD(4,1:2)=COG(b1,:);
                ABCD(2,1:2)=COG(c1,:);
                ABCD(3,1:2)=COG(a1,:);
            elseif (dis1(1)<dis1(2))
                ABCD(4,1:2)=COG(c1,:);
                ABCD(2,1:2)=COG(b1,:);
                ABCD(3,1:2)=COG(a1,:);
            end
        end
    im_static(1:4,1:2,counter)=ABCD;
end
save('im_static','im_static');