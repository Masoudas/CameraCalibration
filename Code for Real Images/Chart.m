% n=5;
% for m=0:36
%     if (m<10)
%         f=imread(strcat('E:\Sequences\Plate 5\Plate-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\Plate 5\Plate-0',num2str(n),'-0',num2str(m),'.jpg'));
%     end
%     f(1:480,623:640,1:3)=1;
%     f(1:31,539:640,1:3)=1;
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\Plate 5\Plate-0',num2str(n),'-00',num2str(m),'.jpg'));
%    elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\Plate 5\Plate-0',num2str(n),'-00',num2str(m),'.jpg'));
%    end
% end

for m=136:340
    f=imread(strcat('E:\Sequences\Plate 6\Plate-06-',num2str(m),'.jpg'));
    g=rgb2gray(f);
    [R C]=size(g);
    level = graythresh(g);
    gb = im2bw(g, 0.5);
    [L,BP] = bwlabel(gb);
    gb=1-gb;
    imwrite(gb,strcat('E:\Sequences\Plate 6\Plate-16-',num2str(m-136),'.jpg'))
end