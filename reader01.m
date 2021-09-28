clc
clear


f=imread('D:\02 The Phd Program in UTDallas\alex\10x.jpg'); % image read
f=im2double(f);  % transform image data into real number
[l w]=size(f(:,:,1)); 
f=rgb2gray(f); % transform image into black and white formate

%diameter=18/2;  
[centers, radii]=imfindcircles(f,[7 10],'ObjectPolarity','bright','Sensitivity',0.98); % find radius and location for each well
r= mean(radii); 
[Nw m]=size(radii);

for i=1:Nw
    rs(i)=r;
    cs(i,1)=round(centers(i,1));
    cs(i,2)=round(centers(i,2));
end 


imshow(f)
h=viscircles(cs, rs)
% intensity=zeros(Nw,1);
%  for k=1:Nw
%      kk=1;
%      for i=1:l
%          for j=1:w
%              if (i-cs(k,2))^2+(j-cs(k,1))^2 < r^2
%                  intensity(k)=intensity(k)+f(i,j); %calculate intensity of the pixel for a well
%                  kk=kk+1;
%              end
%          end
%      end
%      intensity(k)=intensity(k)/kk; % calculate the intensity of a well
%  end
               
                 
     
