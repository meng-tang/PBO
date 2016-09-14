energies =[];
disp('energies for lambdaCurvature from 19 to 30: ');
for lambdaCurvature = 19:30 %19;%0.01%2;%0.1;%1;
lambdaPotts = 0; %0.65
times = [];
lambdaApp = 0.01;
padSize = 3;
bg_mean = 0;
fg_mean = 150;
sigma = 1;
beta = 0; % contrast sensitivity is turned off for now
grayScalar = 0.2;
labels = [1 2];
nType = 4;

%% load the mat file
origImg = double(imread('Stenosis2_Cropped.bmp'));
origImg = double(origImg(:,:,1));
img = origImg;
[numRows, numCols] = size(origImg);

% save(['stenosis2_' num2str(lambdaCurvature) '_results.mat'],'img');

optimizationOptions = setOptimizationOptions_stenosis(grayScalar, ['stenosis2_' num2str(lambdaCurvature)], numRows, numCols, lambdaPotts, lambdaCurvature,lambdaApp, bg_mean, fg_mean, sigma, beta, padSize, nType);


[appUnaryTerms, curvatureUnaryTerms, unaryTerms, subPairwiseTerms, superPairwiseTerms]  = getCurvatureEnergy(numRows, numCols, img, optimizationOptions);
% Meng's wrapper
unaryterms = unaryTerms;
capsource = unaryterms(1,:);
capsink = unaryterms(2,:);
C = sum(capsource);
U = (capsink - capsource)';
pairwiseterms = superPairwiseTerms - subPairwiseTerms;
pair_where = find(pairwiseterms);
pairs = pairwiseterms(pair_where);
pairs_num = size(pairs,1);
M = zeros(pairs_num,3);
img_size = 196*120;
for i=1:pairs_num
    pos = pair_where(i)-1;
    M(i,1:2) = [mod(pos,img_size)+1 int32(floor(pos/img_size)+1)];
    M(i,3) = pairs(i);
end
save('qpdata.mat','C','U','M');
exepath = 'D:/repositories/PBO_public_code/Release/pbo.exe ';
qpdatapath = ' D:/repositories/PBO_public_code/pairwise/qpdata.mat';
[status, result] = dos([exepath qpdatapath ' 196 120']);
%disp(result);
pos=strfind(result,'%');
energy = str2num(result(pos(1)+1:pos(2)-1));
energies = [energies energy];
disp(num2str(energy));
end
return;