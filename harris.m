function harris()
%HARRIS Summary of this function goes here
%   Detailed explanation goes here
k = 0.04;
window_size = 11;
nbpoints = 81;

[file_name, pathname] = uigetfile('*.*','Select path of image file');
fileLocation = strcat(pathname,file_name);
im = imread(fileLocation);
im_orig = im; 
if size(im,3)>1 
    im=rgb2gray(im); 
end 
% Derivative masks 
dx = [-1 0 1; 
      -1 0 1; 
      -1 0 1]; 
dy = dx'; 

% Image derivatives 
Ix = conv2(double(im), dx, 'same'); 
Iy = conv2(double(im), dy, 'same'); 
sigma=2; 
% Generate Gaussian filter of size 9x9 and std. dev. sigma. 
g = fspecial('gaussian',9, sigma); 
% Smoothed squared image derivatives 
Ix2 = conv2(Ix.^2, g, 'same'); 
Iy2 = conv2(Iy.^2, g, 'same'); 
Ixy = conv2(Ix.*Iy, g, 'same');
%figure(1);title();
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1049, 600],'Name',['Corner Detection ', file_name]);
subplot(2,4,1);imshow(im,[]);title(file_name);
hold on

%%Part 1
%Modify the code above to compute a matrix E which contains for every point 
%the value of the smaller eigenvalue of M. 
tic
mask = 1/9*[1 1 1 ;1 1 1 ;1 1 1];
Ix2_c = conv2(Ix2, mask, 'same'); 
Iy2_c = conv2(Iy2, mask, 'same'); 
Ixy_c = conv2(Ixy, mask, 'same');
[i_height, i_width] = size(im);
E = zeros( i_height, i_width);
tic
 for i = 1:i_height
    for j = 1:i_width
        M = [Ix2_c(i,j) Ixy_c(i,j); Ixy_c(i,j) Iy2_c(i,j)];
        eig_v = min(eig(M));
        E(i,j) = eig_v;
    end
 end
elapsed_time_E = toc;
subplot(2,4,3);imshow(mat2gray(E),[]);title(['E -Time= ',num2str(elapsed_time_E),' s']);

%%Part 2

tic
R_loop = zeros(i_height, i_width);
for i = 1:i_height
    for j = 1:i_width
        M = [Ix2_c(i,j), Ixy_c(i,j); Ixy_c(i,j), Iy2_c(i,j)];
        R_loop(i,j) = double(det(M) - k*trace(M)^2);
    end
end
elapsed_time_R_loop = toc;
subplot(2,4,2);imshow(mat2gray(R_loop),[]);title(['R-loop -Time= ',num2str(elapsed_time_R_loop),' s']); 


tic
R = (Ix2_c.* Iy2_c+ Ixy_c.*Ixy_c) - k*(Ix2_c +Iy2_c).^2;
elapsed_time_R_optimized = toc;
subplot(2,4,4);imshow(mat2gray(R),[]);title(['R-optimized -Time= ',num2str(elapsed_time_R_optimized),' s']); 

%%Part 3
tic
features_E_Max = struct('p_x', zeros(nbpoints,1), 'p_y', zeros(nbpoints, 1));
[~, index_E] = sort( E(:), 'descend');
[row_E,col_E] = ind2sub(size(E), index_E);
for i = 1: nbpoints  
    features_E_Max(i).p_y = row_E(i);   
    features_E_Max(i).p_x = col_E(i);  
end
toc
subplot(2,4,5);imshow(im_orig,[]);title('Max - E'); 
hold on;
for i = 1: size(features_E_Max, 2)    
    plot(features_E_Max(i).p_x, features_E_Max(i).p_y, 'r+'); 
end

tic
features_R_Max = struct('p_x', zeros(nbpoints, 1), 'p_y', zeros(nbpoints, 1));
[~, index_R] = sort( R(:), 'descend');
[row_R,col_R] = ind2sub(size(R), index_R);
for i = 1: nbpoints  
    features_R_Max(i).p_y = row_R(i);   
    features_R_Max(i).p_x = col_R(i);  
end
toc

subplot(2,4,6);imshow(im_orig,[]);title('Max - R'); 
hold on;
for i = 1: size(features_R_Max, 2)    
    plot(features_R_Max(i).p_x, features_R_Max(i).p_y, 'r+'); 
end
title('R - Maximal suppresion')


%%Part 4


%Non-maximum supression is often used along with edge 
features_E_NonM = struct('p_x', zeros(nbpoints, 1), 'p_y', zeros(nbpoints, 1));
wsize_2 = floor(window_size/2);
E_pad = padarray(E,[wsize_2,wsize_2]);

count = 1;
while(count<nbpoints)
    for i= 1:size(row_E)  
        if( E_pad(row_E(i),col_E(i)) ~= 0 )
            E_pad(row_E(i)-wsize_2 : row_E(i)+wsize_2, col_E(i)-wsize_2:col_E(i)+wsize_2) = 0;
            features_E_NonM(count).p_y = row_E(i);   
            features_E_NonM(count).p_x = col_E(i);
            count = count + 1;       
        end
        if count == 82
            break;
        end
    end
end

subplot(2,4,7);imshow(im_orig,[]);title('Non Max- E');
hold on;
for i = 1: size(features_E_NonM, 2)    
   plot(features_E_NonM(i).p_x, features_E_NonM(i).p_y, 'r+'); 
   hold on
end
hold on
SubpixAccuracyE.p_x =zeros(nbpoints,1);
SubpixAccuracyE.p_y =zeros(nbpoints,1);
p_x = extractfield(features_E_NonM,'p_x');
p_y = extractfield(features_E_NonM,'p_y');
[SubpixAccuracyE.p_y,SubpixAccuracyE.p_x] = subPixelsAccuracy(E, p_y,p_x,1);


for i=1:nbpoints
plot(SubpixAccuracyE.p_x(i), SubpixAccuracyE.p_y(i), 'g+');
hold on
end


features_R_NonM = struct('p_x', zeros(81, 1), 'p_y', zeros(nbpoints, 1));
R_pad = padarray(R,[wsize_2,wsize_2]);
count = 1;
 while(count<nbpoints)
    for i= 1:size(row_R)  
        if( R_pad(row_R(i),col_R(i)) ~= 0 )
            R_pad(row_R(i)-wsize_2 : row_R(i)+wsize_2, col_R(i)-wsize_2:col_R(i)+wsize_2) = 0;
            features_R_NonM(count).p_y = row_R(i);   
            features_R_NonM(count).p_x = col_R(i);
            count = count + 1;       
        end
        if count == 82
            break;
        end
    end
end
subplot(2,4,8);imshow(im_orig,[]);title('Non Max - R');
hold on;
for i = 1: size(features_R_NonM, 2)    
    plot(features_R_NonM(i).p_x, features_R_NonM(i).p_y, 'r+'); 
end
hold on
SubpixelAccuracyR.p_x =zeros(nbpoints,1);
SubpixelAccuracyR.p_y =zeros(nbpoints,1);
p_x = extractfield(features_R_NonM,'p_x');
p_y = extractfield(features_R_NonM,'p_y');

[SubpixelAccuracyR.p_y,SubpixelAccuracyR.p_x] = subPixelsAccuracy(R, p_y,p_x,1);

hold on
for i=1:nbpoints
plot(SubpixelAccuracyR.p_x(i), SubpixelAccuracyR.p_y(i), 'g+');
end

end

