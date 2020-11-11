clc;
clear all;
close all;

I = double(imread('lena2.tif'));
if ndims(I) > 2
    I = rgb2gray(I);
end

T = 2;
data = randi([0, 1], 1, length(I(:)));
I1 = I; I2 = I;

crossEC1 = zeros(1, 4);
crossEC2 = zeros(1, 4);
crossextractEC1 = zeros(1, 4);
crossextractEC2 = zeros(1, 4);

%-----Embedding in cross: Round 1-----%
[ICrossPred, ec, pc] = crossPredictionDual(I1, 0);
[I1ce, crossEC1(1)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
[ICrossPred, ec, pc] = crossPredictionDual(I2, 1);
[I2co, crossEC2(1)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);

%-----Embedding in dot: Round 1-----%
[ICrossPred, ec, pc] = dotPredictionDual(I2co, 0);
[I1de, crossEC1(2)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
[ICrossPred, ec, pc] = dotPredictionDual(I1ce, 1);
[I2do, crossEC2(2)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);

%-----Embedding in cross: Round 2-----%
[ICrossPred, ec, pc] = crossPredictionDual(I2do, 1);
[I1co, crossEC1(3)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
[ICrossPred, ec, pc] = crossPredictionDual(I1de, 0);
[I2ce, crossEC2(3)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);

%-----Embedding in dot: Round 2-----%
[ICrossPred, ec, pc] = dotPredictionDual(I2ce, 1);
[I1do, crossEC1(4)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);
[ICrossPred, ec, pc] = dotPredictionDual(I1co, 0);
[I2de, crossEC2(4)] = EmbeddingHistogramShifting(ICrossPred, data, T, ec, pc);

%--------------------%
%-----Extraction-----%
%--------------------%

%-----Extraction and Recovery of Dot - Round 2-----%
[ICrossPred, ec, pc] = dotPredictionDual(I1do, 1);%Recovery of I2ce
[Ex2ce, crossextractEC1] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);
[ICrossPred, ec, pc] = dotPredictionDual(I2de, 0);%Recovery of I1co
[Ex1co, crossextractEC2] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);

%-----Extraction and Recovery of Cross - Round 2-----%
[ICrossPred, ec, pc] = crossPredictionDual(Ex2ce, 0);%Recovery of I2ce
[Ex1de, crossextractEC1] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);
[ICrossPred, ec, pc] = crossPredictionDual(Ex1co, 1);%Recovery of I1co
[Ex2do, crossextractEC2] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);

%-----Extraction and Recovery of Dot - Round 1-----%
[ICrossPred, ec, pc] = dotPredictionDual(Ex2do, 1);%Recovery of I2ce
[Ex1ce, crossextractEC1] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);
[ICrossPred, ec, pc] = dotPredictionDual(Ex1de, 0);%Recovery of I1co
[Ex2co, crossextractEC2] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);

%-----Extraction and Recovery of Cross - Round 1-----%
[ICrossPred, ec, pc] = crossPredictionDual(Ex2co, 1);%Recovery of I2ce
[Ex2, crossextractEC1] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);
[ICrossPred, ec, pc] = crossPredictionDual(Ex1ce, 0);%Recovery of I1co
[Ex1, crossextractEC2] = ExtractionHistogramShifting(ICrossPred, ec, T, pc);