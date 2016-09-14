
function newOptimizationOptions = setOptimizationOptions(grayScalar, myString, lambdaPotts, beta, padSize, nType, fg_mu, bg_mu, sigma)


if nargin ==5
    nType = 4;
    fg_mu = 0;
    bg_mu = 1;
    sigma = 0.5;
end
% OPTIMIZATION PARAMS

    newOptimizationOptions.padSize = padSize;
    newOptimizationOptions.neighborhoodBeta = beta;
    newOptimizationOptions.LAMBDA_POTTS = lambdaPotts;
    newOptimizationOptions.NEIGHBORHOOD_TYPE = nType;
    newOptimizationOptions.fg_mu = fg_mu;
    newOptimizationOptions.bg_mu = bg_mu;
    newOptimizationOptions.sigma = sigma;
    newOptimizationOptions.LARGE_CONSTANT = 1000;
    
    % fast trust region params
    newOptimizationOptions.MAX_LAMBDA_GEO = 1e5;
    newOptimizationOptions.LAMBDA_GEO = 2000;%double(0.1);
    newOptimizationOptions.PRECISION_COMPARE_GEO_LAMBDA = 1e-5; % used to compare GEO lambda in parametric maxflow
    newOptimizationOptions.LAMBDA_MULTIPLIER = 2;%5;% used for jumps in backtracking;
    newOptimizationOptions.REDUCTION_RATIO_THRESHOLD = 0.25;%0.25; % used to decide whether to increase or decrease lambda using the multiplier
        
    % output and visualization
    newOptimizationOptions.SHOW_FLAG = false;
    newOptimizationOptions.verbose = 0;
    newOptimizationOptions.textOutput = false;
    newOptimizationOptions.myString = myString;
    newOptimizationOptions.SHOW_RESULT_GRAY_SCALAR = grayScalar;
    
    newOptimizationOptions.correctionFlag = false;
    
   

   
    
end