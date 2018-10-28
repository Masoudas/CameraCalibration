% n=5;
% for m=0:2200
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%     for k=1:3
%         for j=1:65
%             for i=524:640
%                 f(j,i,k)=1;
%             end
%         end
%     end
%     for k=1:3
%         for j=236:261
%             for i=404:462
%                 f(j,i,k)=1;
%             end
%         end
%     end
%     for k=1:3
%         for j=214:257
%             for i=496:612
%                 f(j,i,k)=1;
%             end
%         end
%     end
%     if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=6;
% for m=0:2220
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%     for k=1:3
%         for j=308:314
%             for i=155:189
%                 f(j,i,k)=1;
%             end
%         end
%     end
%     if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=4;
% for m=0:2220
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%        f(143:156,377:413,1:3)=1;
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=4;
% for m=869:1083
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%        f=zeros(480,640,3);
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=3;
% for m=0:2220
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%     f(224:241,209:221,1:3)=1;
% 
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=1;
% for m=0:2220
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%     f(106:111,396:401,1:3)=1;
% 
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=1;
% for m=100:2220
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%     f(1:32,563:640,1:3)=1;
% 
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=1;
% for m=1500:2220
%     if (m<10)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%     else
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
%     f(433:444,393:401,1:3)=1;
% 
%    
%    if (m<10)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-000',num2str(m),'.jpg'));
%     elseif (m<100)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-00',num2str(m),'.jpg'));
%     elseif (m<1000)
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-0',num2str(m),'.jpg'));
%         
%     else
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%     end
% end
% n=2;
% for m=1217:2220
%    
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
%   f(70:76,216:221,1:3)=1;
% 
%    
%   
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% end
% n=2;
% for m=1377:1427
%    
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
%   f(1:480,1:640,1:3)=1;
% 
%    
%   
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% end
%  n=2;
% for m=1851:2220
%    
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% f(77:79,225:227,1:3)=1;
% 
%    
%   
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% end
% n=2;
% for m=1937:2081
%    
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% f(1:480,1:640,1:3)=1;
% 
%    
%   
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% end
% n=5;
% for m=1987:2000
%    
%         f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% f(205:207,52:54,1:3)=1;
% 
%    
%   
%         imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
%    
% end
n=5;
for m=2203:2220
   
        f=imread(strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
   
f(1:480,421:640,1:3)=1;

   
  
        imwrite(f,strcat('E:\Sequences\TFrame-0',num2str(n),'-',num2str(m),'.jpg'));
   
end