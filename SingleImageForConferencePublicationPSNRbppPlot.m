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
    data = randi([0 1], 1, length(I(:)));
    [m, n] = size(I);
    T = 2;
% I = I(:);
%-----Embed in Cross Pixels-----%
[ICrossPred, ec, pc] = crossPrediction(I);
[Ic, crossEC(kk)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
dataEmbeddedCross = data(1:crossEC);

%-----Embed in Dot Pixels-----%
[IDotPred, ed, pd] = dotPrediction(Ic);
[Istego, dotEC(kk)] = EmbeddingHistogramShifting(IDotPred, data(crossEC+1:end), T, ed, pd);
dataEmbeddedDot = data(crossEC+1:crossEC+dotEC);

totalPayload(kk) = crossEC(kk) + dotEC(kk);
    psnr1(kk) = psnr(Istego, I);
    ssim1(kk) = ssim(Istego, I);    
    bpp(kk) = totalPayload(kk)/(m*n);
    
end

ImgNames = {'Barbara', 'Boat', 'Goldhill', 'Lena', 'Mandrill', 'Mountain',...
    'Peppers', 'Washsat', 'Zelda'};
t = table(psnr1', ssim1', totalPayload', bpp', 'RowNames', ImgNames,...
    'VariableNames', {'PSNR', 'SSIM',...
    'payload','bpp'});
display (t)
% filename = 'ConfResults.xlsx';
% writetable(t, filename, 'WriteRowNames', 1, 'Sheet', 1, 'Range', 'A13');