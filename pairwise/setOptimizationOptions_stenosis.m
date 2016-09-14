function newOptimizationOptions = setOptimizationOptions_stenosis...
    (grayScalar, myString, numRows, numCols, lambdaPotts, lambdaCurvature, lambdaApp, bg_mu, fg_mu, sigma, beta, padSize, nType)


% OPTIMIZATION PARAMS

    newOptimizationOptions.numNodes = numRows * numCols;
    newOptimizationOptions.numRows = numRows;
    newOptimizationOptions.numCols = numCols;
    newOptimizationOptions.padSize = padSize;
    frameMask = false(numRows, numCols);
    frameMask([1:padSize (end-padSize+1):end],:) = true;
    frameMask(:,[1:padSize (end-padSize+1):end]) = true;
    newOptimizationOptions.frameMask = frameMask;
    
    newOptimizationOptions.labels = [1 2];

    newOptimizationOptions.neighborhoodBeta = beta;
    newOptimizationOptions.LAMBDA_POTTS = lambdaPotts;
    newOptimizationOptions.LAMBDA_CURVATURE = lambdaCurvature;
    newOptimizationOptions.LAMBDA_APP = lambdaApp;
    newOptimizationOptions.NEIGHBORHOOD_TYPE = nType;
    newOptimizationOptions.fg_mu = fg_mu;
    newOptimizationOptions.bg_mu = bg_mu;
    newOptimizationOptions.sigma = sigma;
    newOptimizationOptions.LARGE_CONSTANT = 1000;
    
    
    % fast trust region params
    newOptimizationOptions.MAX_LAMBDA_GEO = 1e5;
    newOptimizationOptions.LAMBDA_GEO = 2000; %0.1;%double(0.1);
    newOptimizationOptions.PRECISION_COMPARE_GEO_LAMBDA = 1e-5; % used to compare GEO lambda in parametric maxflow
    newOptimizationOptions.LAMBDA_MULTIPLIER = 2;% used for jumps in backtracking;
    newOptimizationOptions.REDUCTION_RATIO_THRESHOLD = 0.25;%0.25; % used to decide whether to increase or decrease lambda using the multiplier
        
    % output and visualization
    newOptimizationOptions.SHOW_FLAG = true;
    newOptimizationOptions.verbose = 0;
    newOptimizationOptions.textOutput = true;
    newOptimizationOptions.myString = myString;
    newOptimizationOptions.SHOW_RESULT_GRAY_SCALAR = grayScalar;
    
    newOptimizationOptions.correctionFlag = false;
    

    
   

   
    
end