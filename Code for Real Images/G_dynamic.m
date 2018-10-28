% Noghte 1, 0 550 0, 2 0 0 0, 3 200 0 0, ...
for m=1:2220
     for n=1
        if (m<10)
            f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
        elseif (m<100)
            f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
        elseif (m<1000)
            f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
        else
            f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
        end
        
        
        g=rgb2gray(f);
        [R C]=size(g);
        level = graythresh(g);
        gb = im2bw(g, 0.6);
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
        if (length(As)<3)
            COGt=-100*(abs(randn(3,2)));
            I=1:3;
        end
        if (length(As)>=4)
            error(msgString)
        end
        for r=1:3
            COG(r,:)=COGt(I(r),:);
        end
        
        % ##############################################################
        % Finding side markers
        % Find A
        D=zeros(2,2);
        dis1(1)=norm(COG(1,:)-COG(2,:));
        dis1(2)=norm(COG(1,:)-COG(3,:));
        dis1(3)=norm(COG(2,:)-COG(3,:));
        [a1, i1]=max(dis1);
        if (i1==1)
            if (dis1(2)>dis1(3))
                D(1,1:2)=COG(1,:);
                D(2,1:2)=COG(2,:);
            elseif (dis1(2)<dis1(3))
                D(1,1:2)=COG(2,:);
                D(2,1:2)=COG(1,:);
            end
        elseif (i1==2)
            if (dis1(1)>dis1(3))
                D(1,1:2)=COG(1,:);
                D(2,1:2)=COG(3,:);
            elseif (dis1(1)<dis1(3))
                D(1,1:2)=COG(3,:);
                D(2,1:2)=COG(1,:);
            end
        elseif (i1==3)
            if (dis1(1)>dis1(2))
                D(1,1:2)=COG(2,:);
                D(2,1:2)=COG(3,:);
            elseif (dis1(1)<dis1(2))
                D(1,1:2)=COG(3,:);
                D(2,1:2)=COG(2,:);
            end
        end
        im(1:2,1:2,m,n)=D;
    end
end

save('im','im');