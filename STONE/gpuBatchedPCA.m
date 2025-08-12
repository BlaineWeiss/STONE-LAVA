function [coeff, explained, meanVec] = gpuBatchedPCA(data, batchSize, numComponents)
% GPUBATCHEDPCA - PCA with GPU acceleration and chunked processing
%
% INPUTS:
%   data          : [N x D] matrix (rows = samples, columns = features)
%   batchSize     : number of rows to process at a time
%   numComponents : number of principal components to return
%
% OUTPUTS:
%   coeff      : [D x numComponents] principal component directions
%   explained  : [numComponents x 1] variance explained (%)
%   meanVec    : [1 x D] mean of the original data

    if nargin < 3
        numComponents = min(50, size(data, 2));
    end

    [N, D] = size(data);

    % Step 1: Compute global mean (CPU)
    fprintf("Computing mean over all data (in batches)...\n");
    meanSum = zeros(1, D);
    totalRows = 0;

    for i = 1:batchSize:N
        idx = i:min(i+batchSize-1, N);
        batch = data(idx, :);
        meanSum = meanSum + sum(batch, 1);
        totalRows = totalRows + size(batch, 1);
    end

    meanVec = meanSum / totalRows;

    % Step 2: Accumulate covariance matrix on GPU
    fprintf("Accumulating covariance matrix on GPU...\n");
    covMat = gpuArray.zeros(D, D);

    for i = 1:batchSize:N
        idx = i:min(i+batchSize-1, N);
        batch = data(idx, :);
        batch_gpu = gpuArray(batch - meanVec);  % center batch
        covMat = covMat + batch_gpu' * batch_gpu;
    end

    covMat = covMat / (totalRows - 1);  % final covariance

    % Step 3: Eigendecomposition (on CPU for stability)
    fprintf("Computing eigendecomposition...\n");
    covMat_cpu = gather(covMat);
    [V, S] = eig(covMat_cpu);

    % Sort eigenvalues & eigenvectors in descending order
    eigVals = diag(S);
    [eigValsSorted, idx] = sort(eigVals, 'descend');
    coeff = V(:, idx(1:numComponents));

    % Explained variance
    totalVar = sum(eigVals);
    explained = 100 * eigValsSorted(1:numComponents) / totalVar;
end