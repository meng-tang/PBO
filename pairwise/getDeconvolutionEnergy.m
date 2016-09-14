% function [unaryTerms pairwiseTerms] = getDeconvolutionEnergy(numRows, numCols, img, optimizationOptions)
% 
%     myFilter = optimizationOptions.myFilter;   
% 
%     filterSize = size(myFilter,1);
%     padSize = filterSize;
%     unaryTerms = zeros(numRows, numCols);
%     myFilterSquared = myFilter(:)*myFilter(:)';
%     myFilterSquared = myFilterSquared(:);
%    
%     pairwiseTerms = sparse(numRows*numCols, numRows*numCols);%, numRows*numCols*(filterSize^4-filterSize^2));
% 
%     
%     counter = 0;
%     for r = (padSize+1):(numRows-padSize)
%         for c = (padSize+1):(numCols-padSize)
%             rowsIdx = (r-floor(filterSize/2)):(r+floor(filterSize/2));
%             colsIdx = (c-floor(filterSize/2)):(c+floor(filterSize/2));
%             unaryTerms(rowsIdx,colsIdx) = ...
%                 unaryTerms(rowsIdx,colsIdx)...
%                 -2*img(r,c).*myFilter;
%                 
%             [frows fcols] = ndgrid((r-floor(filterSize/2)):(r+floor(filterSize/2)),(c-floor(filterSize/2)):(c+floor(filterSize/2)));
%             nodeIds = frows + (fcols-1)*numRows;
%             nodeIds = nodeIds(:);
%             [nodeIds neighborNodeIds] = ndgrid(nodeIds,nodeIds);
%             
%             pairwiseTerms = pairwiseTerms + sparse( [nodeIds(:); neighborNodeIds(:)],[neighborNodeIds(:); nodeIds(:)], [myFilterSquared/2; myFilterSquared/2], numRows*numCols, numRows*numCols);
%               
%             counter = counter + 1;
%             if ~mod(r*c,500)
%                 disp(num2str(counter/(numRows*numCols)));
%             end
%         end
%     end
%     
%     
%     unaryTerms = unaryTerms(:)' + pairwiseTerms(logical(speye(numRows*numCols)))';
%     unaryTerms = [zeros(1,numRows*numCols); unaryTerms];
%     pairwiseTerms(logical(speye(numRows*numCols))) = 0;
% %     pairwiseTerms = sparse(tempNodeIds, tempNeighborIds, tempWeights,numRows*numCols, numRows*numCols, numRows*numCols*(filterSize^4-filterSize^2));
%     
%     
% end


function [unaryTerms, pairwiseTerms, constantTerm] = getDeconvolutionEnergy(numRows, numCols, img, optimizationOptions)

    myFilter = optimizationOptions.myFilter;  
    filterSize = size(myFilter,1);
    center = ceil(filterSize/2);
        
    constantTerm = 0;
    
    % unary terms
    unaryTerms = zeros(numRows, numCols);
    for r = center:(numRows-center+1)
        for c = center:(numCols-center+1)
            rowsIdx = (r-floor(filterSize/2)):(r+floor(filterSize/2));
            colsIdx = (c-floor(filterSize/2)):(c+floor(filterSize/2));
            unaryTerms(rowsIdx,colsIdx) = ...
                unaryTerms(rowsIdx,colsIdx)...
                -2*img(r,c).*myFilter;
            constantTerm = constantTerm + img(r,c)^2;
        end
    end
    
    
    % supermodular pairwise terms
    pairwiseTerms = sparse(numRows*numCols, numRows*numCols);    
    [shiftRows, shiftCols] = ndgrid((1:filterSize)-center,(1:filterSize)-center);
    
   
    [allRows, allCols] = ndgrid(1:numRows, 1:numCols);
    % this idxImg will hold site id for each pixel (rowwise)
    idxImg = zeros(numRows, numCols);
    
    % idx scans img colwise
    idxImg(:) = 1:(numRows*numCols);
    
    for pShift = 1:numel(shiftRows)
        pRowShift = shiftRows(pShift);
        pColShift = shiftCols(pShift);
        pWeight = myFilter(pShift);
    
        for qShift = 1:numel(shiftRows)
            qRowShift = shiftRows(qShift);
            qColShift = shiftCols(qShift);
            qWeight = myFilter(qShift);
    
            % add the shift and only take those pixels that have neighbors
            % with this shift inside the image
            pRows = allRows + pRowShift;
            idx = (pRows <= numRows) & (pRows>=1);

            pCols = allCols + pColShift;
            idx = (pCols <= numCols) & (pCols>=1) & idx;


            qRows = allRows + qRowShift;
            idx = (qRows <= numRows) & (qRows>=1) & idx;

            qCols = allCols + qColShift;
            idx = (qCols <= numCols) & (qCols>=1) & idx;


            % these are the indices of the pairs of pixels 
            %p = vec(idxImg(sub2ind([numRows, numCols],pRows(idx),pCols(idx))));
            %q = vec(idxImg(sub2ind([numRows, numCols],qRows(idx),qCols(idx))));
            
            p = idxImg(sub2ind([numRows, numCols],pRows(idx),pCols(idx)));
            q = idxImg(sub2ind([numRows, numCols],qRows(idx),qCols(idx)));

            % make it summetric pairwise potential
            siteIdx = [p q];
            nSiteIdx = [q p];
            weights = [0.5*pWeight*qWeight*ones(length(p),1) 0.5*pWeight*qWeight*ones(length(p),1)];

            pairwiseTerms = pairwiseTerms + sparse( siteIdx, nSiteIdx, weights, numRows*numCols, numRows*numCols);    
        end
    end

    % diagonal pairwise terms are just linear terms - remove them and add
    % to unary
    unaryTerms = unaryTerms(:)' + pairwiseTerms(logical(speye(numRows*numCols)))';
    pairwiseTerms(logical(speye(numRows*numCols))) = 0;
    unaryTerms = [zeros(1,numRows*numCols); unaryTerms];
    
end