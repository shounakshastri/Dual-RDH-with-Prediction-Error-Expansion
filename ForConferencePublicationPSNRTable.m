clc;
clear all;
close all;

imgSet=imageSet('C:\Users\shoun\Documents\MATLAB\Image Databases\Waterloo');
no=imgSet.Count;
psnr1=zeros(1,no);
psnr2=zeros(1,no);
ssim1=zeros(1,no);
ssim2=zeros(1,no);
totalPayload=zeros(1,no);
bpp=zeros(1,no);

for kk = 1:no
    kk
    I=double(imgSet.read(kk));
    [m, n] = size(I);
data = randi([0, 1], 1, m*n);
l = length(data);
I = padarray(I, [1,1]);
p1=zeros(size(I)); p2=zeros(size(I));
I1 = I; I2 = I;
T = 2;

% Cross Embedding Round 1%
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii, 2) == 0
            if mod(ii+jj, 2) == 0
                I1(ii, jj) = floor((I1(ii+1, jj) + I1(ii-1, jj) + I1(ii, jj+1) + I1(ii, jj-1))/4);
                p1(ii,jj)=1;
            end
        else
            if mod(ii+jj, 2) == 0
                I2(ii, jj) = floor((I2(ii+1, jj) + I2(ii-1, jj) + I2(ii, jj+1) + I2(ii, jj-1))/4);
                p2(ii,jj)=1;
            end
        end
    end
end
I = I(2:end-1, 2:end-1);
I1 = I1(2:end-1, 2:end-1); I2 = I2(2:end-1, 2:end-1);
p1 = p1(2:end-1, 2:end-1); p2 = p2(2:end-1, 2:end-1);
er1 = I - I1; er2 = I - I2;

ed1 = zeros(m, n); ed2 = zeros(m, n);%Shifted Prediction Errors
temp1 = 0; temp2 = 0;
for ii = 1:m
    for jj = 1:n
        if l <= temp1 + temp2
            ed1(ii, jj) = er1(ii, jj);
            ed2(ii, jj) = er2(ii, jj);
        else
            if p1(ii, jj) == 1
                if er1(ii, jj) >= -T && er1(ii, jj) < T
                    temp1 = temp1 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed1(ii, jj) = 2*er1(ii, jj) + data(temp1 + temp2);                
                elseif er1(ii, jj) >= T
                    ed1(ii, jj) = er1(ii, jj) + T;
                elseif er1(ii, jj) < -T
                    ed1(ii, jj) = er1(ii, jj) - T;
                end
            elseif p2(ii, jj) == 1
                if er2(ii, jj) >= -T && er2(ii, jj) < T
                    temp2 = temp2 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed2(ii, jj) = 2*er2(ii, jj) + data(temp1 + temp2);                
                elseif er2(ii, jj) >= T
                    ed2(ii, jj) = er2(ii, jj) + T;
                elseif er2(ii, jj) < -T
                    ed2(ii, jj) = er2(ii, jj) - T;
                end
            end
        end
    end
end
I_stego1 = I1 + ed1; I_stego2 = I2 + ed2;

% Dot Embedding Round 1%
I1 = I_stego1; I2 = I_stego2;
I1 = padarray(I1, [1, 1]); I2 = padarray(I2, [1, 1]);
p1 = zeros(size(I1)); p2 = zeros(size(I2));
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii, 2) == 0
            if mod(ii+jj, 2) ~= 0
                I1(ii, jj) = floor((I1(ii+1, jj) + I1(ii-1, jj) + I1(ii, jj+1) + I1(ii, jj-1))/4);
                p1(ii,jj)=1;
            end
        else
            if mod(ii+jj, 2) ~= 0
                I2(ii, jj) = floor((I2(ii+1, jj) + I2(ii-1, jj) + I2(ii, jj+1) + I2(ii, jj-1))/4);
                p2(ii,jj)=1;
            end
        end
    end
end

I1 = I1(2:end-1, 2:end-1); I2 = I2(2:end-1, 2:end-1);
p1 = p1(2:end-1, 2:end-1); p2 = p2(2:end-1, 2:end-1);
er1 = I - I1; er2 = I - I2;

ed1 = zeros(m, n); ed2 = zeros(m, n);%Shifted Prediction Errors

for ii = 1:m
    for jj = 1:n
        if l <= temp1 + temp2
            ed1(ii, jj) = er1(ii, jj);
            ed2(ii, jj) = er2(ii, jj);
        else
            if p1(ii, jj) == 1
                if er1(ii, jj) >= -T && er1(ii, jj) < T
                    temp1 = temp1 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed1(ii, jj) = 2*er1(ii, jj) + data(temp1 + temp2);                
                elseif er1(ii, jj) >= T
                    ed1(ii, jj) = er1(ii, jj) + T;
                elseif er1(ii, jj) < -T
                    ed1(ii, jj) = er1(ii, jj) - T;
                end
            elseif p2(ii, jj) == 1
                if er2(ii, jj) >= -T && er2(ii, jj) < T
                    temp2 = temp2 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed2(ii, jj) = 2*er2(ii, jj) + data(temp1 + temp2);                
                elseif er2(ii, jj) >= T
                    ed2(ii, jj) = er2(ii, jj) + T;
                elseif er2(ii, jj) < -T
                    ed2(ii, jj) = er2(ii, jj) - T;
                end
            end
        end
    end
end
I_stego1 = I1 + ed1; I_stego2 = I2 + ed2;

% Cross Embedding Round 2

I1 = I_stego1; I2 = I_stego2;
I1 = padarray(I1, [1, 1]); I2 = padarray(I2, [1, 1]);
p1 = zeros(size(I1)); p2 = zeros(size(I2));
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii, 2) ~= 0
            if mod(ii+jj, 2) == 0
                I1(ii, jj) = floor((I1(ii+1, jj) + I1(ii-1, jj) + I1(ii, jj+1) + I1(ii, jj-1))/4);
                p1(ii,jj)=1;
            end
        else
            if mod(ii+jj, 2) == 0
                I2(ii, jj) = floor((I2(ii+1, jj) + I2(ii-1, jj) + I2(ii, jj+1) + I2(ii, jj-1))/4);
                p2(ii,jj)=1;
            end
        end
    end
end

I1 = I1(2:end-1, 2:end-1); I2 = I2(2:end-1, 2:end-1);
p1 = p1(2:end-1, 2:end-1); p2 = p2(2:end-1, 2:end-1);
er1 = I - I1; er2 = I - I2;

ed1 = zeros(m, n); ed2 = zeros(m, n);%Shifted Prediction Errors

for ii = 1:m
    for jj = 1:n
        if l <= temp1 + temp2
            ed1(ii, jj) = er1(ii, jj);
            ed2(ii, jj) = er2(ii, jj);
        else
            if p1(ii, jj) == 1
                if er1(ii, jj) >= -T && er1(ii, jj) < T
                    temp1 = temp1 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed1(ii, jj) = 2*er1(ii, jj) + data(temp1 + temp2);                
                elseif er1(ii, jj) >= T
                    ed1(ii, jj) = er1(ii, jj) + T;
                elseif er1(ii, jj) < -T
                    ed1(ii, jj) = er1(ii, jj) - T;
                end
            elseif p2(ii, jj) == 1
                if er2(ii, jj) >= -T && er2(ii, jj) < T
                    temp2 = temp2 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed2(ii, jj) = 2*er2(ii, jj) + data(temp1 + temp2);                
                elseif er2(ii, jj) >= T
                    ed2(ii, jj) = er2(ii, jj) + T;
                elseif er2(ii, jj) < -T
                    ed2(ii, jj) = er2(ii, jj) - T;
                end
            end
        end
    end
end
I_stego1 = I1 + ed1; I_stego2 = I2 + ed2;


% Dot Embedding Round 2 %
I1 = I_stego1; I2 = I_stego2;
I1 = padarray(I1, [1, 1]); I2 = padarray(I2, [1, 1]);
p1 = zeros(size(I1)); p2 = zeros(size(I2));
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii, 2) ~= 0
            if mod(ii+jj, 2) ~= 0
                I1(ii, jj) = floor((I1(ii+1, jj) + I1(ii-1, jj) + I1(ii, jj+1) + I1(ii, jj-1))/4);
                p1(ii,jj)=1;
            end
        else
            if mod(ii+jj, 2) ~= 0
                I2(ii, jj) = floor((I2(ii+1, jj) + I2(ii-1, jj) + I2(ii, jj+1) + I2(ii, jj-1))/4);
                p2(ii,jj)=1;
            end
        end
    end
end

I1 = I1(2:end-1, 2:end-1); I2 = I2(2:end-1, 2:end-1);
p1 = p1(2:end-1, 2:end-1); p2 = p2(2:end-1, 2:end-1);
er1 = I - I1; er2 = I - I2;

ed1 = zeros(m, n); ed2 = zeros(m, n);%Shifted Prediction Errors

for ii = 1:m
    for jj = 1:n
        if l <= temp1 + temp2
            ed1(ii, jj) = er1(ii, jj);
            ed2(ii, jj) = er2(ii, jj);
        else
            if p1(ii, jj) == 1
                if er1(ii, jj) >= -T && er1(ii, jj) < T
                    temp1 = temp1 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed1(ii, jj) = 2*er1(ii, jj) + data(temp1 + temp2);                
                elseif er1(ii, jj) >= T
                    ed1(ii, jj) = er1(ii, jj) + T;
                elseif er1(ii, jj) < -T
                    ed1(ii, jj) = er1(ii, jj) - T;
                end
            elseif p2(ii, jj) == 1
                if er2(ii, jj) >= -T && er2(ii, jj) < T
                    temp2 = temp2 + 1;
                    %                 data(temp) = randi([0, 1], 1);
                    ed2(ii, jj) = 2*er2(ii, jj) + data(temp1 + temp2);                
                elseif er2(ii, jj) >= T
                    ed2(ii, jj) = er2(ii, jj) + T;
                elseif er2(ii, jj) < -T
                    ed2(ii, jj) = er2(ii, jj) - T;
                end
            end
        end
    end
end
I_stego1 = I1 + ed1; I_stego2 = I2 + ed2;
% Calculations %
    psnr1(kk) = psnr(I_stego1, I);
    psnr2(kk) = psnr(I_stego2, I);
    ssim1(kk) = ssim(I_stego1, I);
    ssim2(kk) = ssim(I_stego2, I);
    totalPayload(kk) = temp1 + temp2;
    bpp(kk) = totalPayload(kk)/(m*n);
    
end

ImgNames = {'Barbara', 'Boat', 'Goldhill', 'Lena', 'Mandrill', 'Mountain',...
    'Peppers', 'Washsat', 'Zelda'};
t = table(psnr1', psnr2', ssim1', ssim2', totalPayload', bpp', 'RowNames', ImgNames,...
    'VariableNames', {'PSNR1', 'PSNR2', 'SSIM1',...
    'SSIM2', 'Payload', 'bpp'});
filename = 'ConfResults.xlsx';
writetable(t, filename, 'WriteRowNames', 1);