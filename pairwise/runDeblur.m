noiseSigma = 0.2;
energies=[];
for k=1:10
lambdaPotts = 0;%0.01; %0.65
filterSize = 3;
padSize = filterSize;
myFilter = ones(filterSize)/filterSize^2;
beta = 0; % contrast sensitivity is turned off for now
grayScalar = 0.2;
labels = [1 2];
myString = ['ECCV_' num2str(noiseSigma)];
optimizationOptions = setOptimizationOptions(grayScalar, myString, lambdaPotts, beta, padSize);
optimizationOptions.myFilter = myFilter;

%% load the mat file
origImg = double(imread('ECCV.png'))/255;
origImg = origImg(:,:,1);
origImg(origImg>0.5) = 1;
origImg(origImg<0.5) = 0;
[numRows, numCols] = size(origImg);
img = imfilter(origImg,myFilter) + normrnd(0,noiseSigma,numRows, numCols);
resultsFileName = [myString '_blur.png'];
imwrite(img,resultsFileName);
% load('CVPR_blur_noise'); % the same random initialization for all the methods to compare
initLabeling = 2*ones(numRows, numCols);
% initLabeling([1:3 end-3:end],[1:3 end-3:end]) = 1;
% initLabeling(floor(numRows/4):(numRows- floor(numRows/4)),floor(numCols/4):(numCols- floor(numCols/4))) = 2;
%%
[unaryTerms, superPairwiseTerms, constantTerm] = getDeconvolutionEnergy(numRows, numCols, img, optimizationOptions);
unaryterms = unaryTerms;
capsource = unaryterms(1,:);
capsink = unaryterms(2,:);
C = sum(capsource)+constantTerm;
U = (capsink - capsource)';
pairwiseterms = superPairwiseTerms;
pair_where = find(pairwiseterms);
pairs = pairwiseterms(pair_where);
pairs_num = size(pairs,1);
M = zeros(pairs_num,3);
img_size = 120*104;
for i=1:pairs_num
    pos = pair_where(i)-1;
    M(i,1:2) = [mod(pos,img_size)+1 int32(floor(pos/img_size)+1)];
    M(i,3) = pairs(i);
end
save('qpdata.mat','C','U','M');
exepath = 'D:/repositories/PBO_public_code/Release/pbo.exe ';
qpdatapath = ' D:/repositories/PBO_public_code/pairwise/qpdata.mat';
[status, result] = dos([exepath qpdatapath ' 120 104']);
%disp(result);
pos=strfind(result,'%');
energy = str2num(result(pos(1)+1:pos(2)-1));
energies = [energies energy];
disp(['run ' num2str(k) ' st time: energy = ' num2str(energy)]);
end
disp(['mean energy: ' num2str(mean(energies))]);
return;