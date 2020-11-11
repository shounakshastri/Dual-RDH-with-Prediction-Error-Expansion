clc;
clear all;
close all;
% I = [46 44 44 38;...
%      49 50 43 34;...
%      53 39 38 32;...
%      45 36 45 46];
I = double(imread('lena2.tif'));
% figure, imshow(I, []);
[m, n] = size(I);

I1 = padarray(I, [1,1]);%This is a padded image. No operations
                                     %are done on this image. DO NOT USE
                                     %THE 'REPLICATE' OPTION HERE !!!
p = zeros(size(I1));%This is the map which indicates the cross pixels by 1.
I2 = I1;%This is the Cross Predicted Image. 
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii+jj, 2) == 0
            I2(ii, jj) = floor((I1(ii+1, jj) + I1(ii-1, jj) + I1(ii, jj+1) + I1(ii, jj-1))/4);
            p(ii,jj)=1;
        end
    end
end
I2 = I2(2:end-1, 2:end-1);%Predicted Image
p = p(2:end-1, 2:end-1);%Map of cross pixels
e_c = I - I2;%Prediction Errors
% figure, imshow(I2, []);

%-----Cross Embedding-----%
e_dash_c = zeros(m, n);%Shifted Prediction Errors
temp = 1;
T = 2;
for ii = 1:m
    for jj = 1:n
        if p(ii, jj) == 1
            if e_c(ii, jj) >= -T && e_c(ii, jj) < T
                b(temp) = randi([0, 1], 1);
                e_dash_c(ii, jj) = 2*e_c(ii, jj) + b(temp);
                temp = temp + 1;
            elseif e_c(ii, jj) >= T
                e_dash_c(ii, jj) = e_c(ii, jj) + T;
            elseif e_c(ii, jj) < -T
                e_dash_c(ii, jj) = e_c(ii, jj) - T;
            end
        end
    end
end
I3 = I2 + e_dash_c;%Stego Image
% figure, imshow(I3, []);

%-----Cross Extraction-----%

I4 = padarray(I3, [1,1]);%This is a padded image. No operations
                                     %are done on this image.

I5 = I4;%This is the Cross Predicted Image. 
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii+jj, 2) == 0
            I5(ii, jj) = floor((I4(ii+1, jj) + I4(ii-1, jj) + I4(ii, jj+1) + I4(ii, jj-1))/4);
        end
    end
end
I6 = I5(2:end-1, 2:end-1);
ex_c = I3 - I6;

temp = 1;
e_c_1 = zeros(m, n);
for ii = 1:m
    for jj = 1:n
        if p(ii, jj) == 1
            if ex_c(ii, jj) >= -2*T && ex_c(ii, jj) < 2*T
                b_1(temp) = mod(ex_c(ii, jj), 2);
                e_c_1(ii, jj) = floor(ex_c(ii, jj)/2);
                temp = temp + 1;
            elseif ex_c(ii, jj) < -2*T 
                e_c_1(ii, jj) = ex_c(ii, jj) + T;
            elseif ex_c(ii, jj) >= 2*T
                e_c_1(ii, jj) = ex_c(ii, jj) - T;
            end
        end
    end
end
I7 = I6 + e_c_1;%Recovered Image
disp([isequal(e_c, e_c_1) isequal(b, b_1) isequal(I, I7)]);
% figure, imshow(I7, []);
