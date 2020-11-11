clc;
clear all;
close all;

% This code is created to test the dual embedding and extraction functions.
% The output is tested using the isequal() function. If both outputs are 1
% (True) then the code is working well.

I = double(imread('lena2.tif'));
T = 2;
data = randi([0, 1], 1, length(I(:)));

%-----Embedding in cross and dot: even rows-----%
[ICrossPred, ec, pc] = crossPredictionDual(I, 0);
[Ic, crossEC1] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
[ICrossPred, ec, pc] = dotPredictionDual(Ic, 0);
[Id, crossEC2] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);

%-----Extraction from cross and dot: even rows-----%
[IdrossPredExtract, ecExtract, pc] = dotPredictionDual(Id, 0);
[Irec, dataCrossRec] = ExtractionHistogramShifting(IdrossPredExtract, ecExtract, T, pc);
[ICrossPredExtract, ecExtract, pc] = crossPredictionDual(Irec, 0);
[Idrec, dataCrossRec] = ExtractionHistogramShifting(ICrossPredExtract, ecExtract, T, pc);

[isequal(I,Idrec) isequal(data(1:crossEC1), dataCrossRec)]
psnr(I, Ic)

%-----Embedding in cross and dot: odd rows-----%
[ICrossPred, ec, pc] = crossPredictionDual(I, 1);
[Ic, crossEC1] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
[ICrossPred, ec, pc] = dotPredictionDual(Ic, 1);
[Id, crossEC2] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);

%-----Extraction cross and dot: odd rows-----%
[IdrossPredExtract, ecExtract, pc] = dotPredictionDual(Id, 1);
[Irec, dataCrossRec] = ExtractionHistogramShifting(IdrossPredExtract, ecExtract, T, pc);
[ICrossPredExtract, ecExtract, pc] = crossPredictionDual(Irec, 1);
[Idrec, dataCrossRec] = ExtractionHistogramShifting(ICrossPredExtract, ecExtract, T, pc);

[isequal(I,Idrec) isequal(data(1:crossEC1), dataCrossRec)]
psnr(I, Ic)
