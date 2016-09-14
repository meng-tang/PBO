function [appUnaryTerms, curvatureUnaryTerms, unaryTerms, subPairwiseTerms, superPairwiseTerms] =...
    getCurvatureEnergy(numRows, numCols, img, optimizationOptions)

nhMaskX1 = [0,1,0,-1]; 
nhMaskY1 = [-1,0,1,0];

nhMaskX3 = [1,0,-1,0]; 
nhMaskY3 = [0,1,0,-1];


nhNum = length(nhMaskX1);
nhWeights = optimizationOptions.LAMBDA_CURVATURE * pi^2/4 * ones(size(nhMaskX1));
    % energy terms look as follows (for each direction alpha)
    % 
    % unary:
    % w*x_2 - 0.5*w*x_1 - 0.5*w*x_2 - 0.5*w*x_2 - 0.5*w*x_3
    %         \-- unary from pw1 --/ \-- unary from pw2 --/
    %    
    % submodular pairwise:
    % + 0.5*w*x_1(1-x_2) + 0.5*w*x_2*(1-x_1)
    % + 0.5*w*x_2(1-x_3) + 0.5*w*x_3*(1-x_2)
    %
    % non-submodular pairwise:\
    % + w*x_1*x_3
    %
    % here x_1, x_2 and x_3 are defined as follows for each direction
    % ----------------
    % |    |    |    |
    % |    |x_1 |    |
    % ----------------
    % |    | |  |    |
    % |    |x_2-|x_3 |
    % ----------------
    % |    |    |    |
    % |    |    |    |
    % ----------------

    % define unary terms which will hold both appearance and curvature
    % unary terms
    appUnaryTerms = zeros(2, numRows * numCols);
    curvatureUnaryTerms  = zeros(2, numRows * numCols);
    
    appUnaryTerms(1,:) = optimizationOptions.LAMBDA_APP*(optimizationOptions.bg_mu - img(:)).^2 ;
    appUnaryTerms(2,:) = optimizationOptions.LAMBDA_APP*(optimizationOptions.fg_mu - img(:)).^2;
    % frameMask
%     appUnaryTerms(2,optimizationOptions.frameMask) = 1000;
%     appUnaryTerms(1,optimizationOptions.frameMask) = 0;
%     
    subPairwiseTerms = sparse(numRows*numCols, numRows*numCols);
    superPairwiseTerms = sparse(numRows*numCols, numRows*numCols);
    
    numPixels = numRows*numCols;
    
    % for all constraints construct unary terms
    for i=1:nhNum
        
            % compute indices for x_1, x_2 and x_3
            [x2Y, x2X] = meshgrid(1:numRows, 1:numCols);

            x1X = x2X + nhMaskX1(i);
            x1Y = x2Y + nhMaskY1(i);

            x3X = x2X + nhMaskX3(i);
            x3Y = x2Y + nhMaskY3(i);
            
            % indices for pairwise terms x_1,x_2: set only those terms inside image
            idx = find(x3X>0 & x3X<=numCols & x3Y>0 & x3Y<=numRows & ...
                       x1X>0 & x1X<=numCols & x1Y>0 & x1Y<=numRows);

            idxX1 = sub2ind([numRows, numCols], x1Y(idx), x1X(idx));
            % use only those submodular pairwise terms where x_1 is inside image
            idxX2 = sub2ind([numRows, numCols], x2Y(idx), x2X(idx)); %indices for x2 for pairwise between 1 and 2
            % indices for pairwise terms x_2,x_3: set only those terms inside image
            idxX3 = sub2ind([numRows, numCols], x3Y(idx), x3X(idx));

            
            % set supermodular pairwise term between x_1 and x_3
            % we are making it symmertric: both (x1, x3) and (x3,x1), that is why 
            % we multiply by 0.5
            superPairwiseTerms = superPairwiseTerms + sparse( ...
                [idxX3 idxX1], [idxX1 idxX3], 0.5* nhWeights(i), numPixels, numPixels);
            

            % set the submodular terms
            % (x_1, x_2)
            % set unary terms
            %curvatureUnaryTerms(2,idxX1) = curvatureUnaryTerms(2,idxX1) - 0.5*nhWeights(i);
            % potts pairwise term between x_1 and x_2 (symmetric)
            %subPairwiseTerms = subPairwiseTerms + sparse( ...
                %[idxX1, idxX2], [idxX2 idxX1], 0.5*nhWeights(i), numPixels, numPixels);
            % (x_2, x_3)
            % set unary terms
            %curvatureUnaryTerms(2,idxX3) = curvatureUnaryTerms(2,idxX3) - 0.5*nhWeights(i);
            
            % potts pairwise term between x_2 and x_3
            %subPairwiseTerms = subPairwiseTerms + sparse( ...
                %[idxX3 idxX2], [idxX2 idxX3], 0.5*nhWeights(i), numPixels, numPixels);
            
            % Meng
            curvatureUnaryTerms(2,idxX2) = curvatureUnaryTerms(2,idxX2) + nhWeights(i);
            subPairwiseTerms = subPairwiseTerms + sparse( ...
                [idxX3 idxX2], [idxX2 idxX3], 0.5*nhWeights(i), numPixels, numPixels);
            subPairwiseTerms = subPairwiseTerms + sparse( ...
                [idxX1 idxX2], [idxX2 idxX1], 0.5*nhWeights(i), numPixels, numPixels);
            
            
            
    end  
    unaryTerms = curvatureUnaryTerms + appUnaryTerms;
    
end