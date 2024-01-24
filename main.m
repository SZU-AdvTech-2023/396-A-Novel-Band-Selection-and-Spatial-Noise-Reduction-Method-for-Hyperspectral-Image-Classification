%This repository contains code developed for the Computer Frontiers Technology course assignment at Shenzhen University. The assignment was completed by student Zhang Yuze, with student ID 2350278003. The code is a reproduction of the research paper titled "A Novel Band Selection and Spatial Noise Reduction Method for Hyperspectral Image Classification."


groundTruth = allLabels;
imageData = allData;
[Nx, Ny, numBands] = size(imageData);

%% Eliminate zero values
reshapedData = reshape(imageData, Nx * Ny, numBands);
reshapedData = reshapedData(:);
zeroIndices = find(reshapedData == 0);
reshapedData(zeroIndices) = 0.001;
normalizedData = reshape(reshapedData, Nx * Ny, numBands);
imageData = reshape(normalizedData, Nx, Ny, numBands);

tic;

%% NGNMF for band selection
reshapedData = reshape(imageData, Nx * Ny, numBands);
reshapedData = reshapedData';
reshapedData = reshapedData(:);
reshapedData = mapminmax(reshapedData', 0, 1);  % normalization
reshapedData(find(reshapedData) == 0) = 0.001;
normalizedData = reshape(reshapedData, numBands, Nx * Ny);
similarityMatrix = computeSimilarityMatrix(normalizedData');  % compute similarity matrix
selectedBands = N_band;  % selected band number
K = selectedBands + 1;
bandSubspace = computeSelectedBands(similarityMatrix, numBands, K);

% Band Ranking With the Normalization MF
for i = 1:numBands
    normalizedImage(:,:,i) = mapminmax(imageData(:,:,i), 0, 1);
end
reshapedNormalizedImage = reshape(normalizedImage, Nx * Ny, numBands);
[bandOrder] = NMFW_RANK(reshapedNormalizedImage);
resultBandId = zeros(1, selectedBands);
for i = 1:selectedBands
    [~, ind] = ismember(bandSubspace(i):bandSubspace(i + 1), bandOrder);
    resultBandId(i) = find(ind == min(ind)) + bandSubspace(i) - 1;
end
selectedImage = imageData(:,:,resultBandId);

%% E2DSSA for Spatial FE
u = 3; w = 2 * u + 1; w2 = w * w;  % search window
L = 25;  % number of embedded pixels
for i = 1:numBands
    paddedImage(:,:,i) = padarray(imageData(:,:,i), [u, u], 'symmetric', 'both');
end
moderesultAN = zeros(Nx * w, Ny * w);
for i = 1:Nx
    for j = 1:Ny
        i1 = i + u; j1 = j + u;
        testCube = paddedImage(i1 - u:i1 + u, j1 - u:j1 + u, :);
        m = reshape(testCube, [w2, numBands]);
        center = m((w2 + 1) / 2, :);
        NED = zeros(1, w2);
        for ii = 1:w2
            NED(:, ii) = sqrt(sum(power((m(ii, :) / norm(m(ii, :)) - center / norm(center)), 2)));  % NED
        end
        [~, ind] = sort(NED);
        index = ind(1:L);

        mask = zeros(w2, 1);
        mask(index) = 1;
        mask = reshape(mask, w, w);
        moderesultAN(w * (i - 1) + 1:w * i, w * (j - 1) + 1:w * j) = mask;
    end
end

% Band-Based Spatial Processing
spatialImage = zeros(Nx, Ny, selectedBands);
selectedImagePadding = paddedImage(:,:,resultBandId);
for k = 1:selectedBands
    X2D = zeros(L, Nx * Ny);
    ID = zeros(L, Nx * Ny);
    % Adaptive embedding
    n = 0;
    for i = 1:Nx
        for j = 1:Ny
            n = n + 1;
            i1 = i + u; j1 = j + u;
            test = selectedImagePadding(i1 - u:i1 + u, j1 - u:j1 + u, k);
            col = moderesultAN(w * (i - 1) + 1:w * i, w * (j - 1) + 1:w * j) .* test;
            col = col(:);
            selCol = col(col ~= 0);
            indexCol = find(col ~= 0);
            X2D(:, n) = selCol;
            ID(:, n) = indexCol;
        end
    end
    % Singular Value Decomposition and grouping
    S = X2D * X2D';
    [U, autoval] = eigs(S, 1);
    V = (X2D') * U;
    rca = U * V';
    % Reprojection
    newPaddedImage = zeros(Nx + w - 1, Ny + w - 1);
    repeat = zeros(Nx + w - 1, Ny + w - 1);
    kk = 0;
    for i = 1:Nx
        for j = 1:Ny
            kk = kk + 1;
            recCol = zeros(w2, 1);
            recCol(ID(:, kk)) = rca(:, kk);

            i1 = i + u; j1 = j + u;
            newPaddedImage(i1 - u:i1 + u, j1 - u:j1 + u) = newPaddedImage(i1 - u:i1 + u, j1 - u:j1 + u) + reshape(recCol, w, w);
            repeat(i1 - u:i1 + u, j1 - u:j1 + u) = repeat(i1 - u:i1 + u, j1 - u:j1 + u) + moderesultAN(w * (i - 1) + 1:w * i, w * (j - 1) + 1:w * j);
        end
    end
    newPaddedImage = newPaddedImage ./ repeat;
    spatialImage(:,:,k) = newPaddedImage(u + 1:Nx + u, u + 1:Ny + u);
end
toc;

%% Training-test samples
labels = groundTruth(:);
vectors = reshape(spatialImage, Nx * Ny, selectedBands);
classNum = max(max(groundTruth)) - min(min(groundTruth));
trainVectors = [];
trainLabels = [];
trainIndex = [];
testVectors = [];
testLabels = [];
testIndex = [];
rng('default');
sampPro = 0.8;  % proportion of training samples
for k = 1:1:classNum
    index = find(labels == k);
    perclassNum = length(index);
    vectorsPerclass = vectors(index, :);
    c = randperm(perclassNum);
    selectTrain = vectorsPerclass(c(1:ceil(perclassNum * sampPro)), :);  % select training samples
    trainIndexK = index(c(1:ceil(perclassNum * sampPro)));
    trainIndex = [trainIndex; trainIndexK];
    selectTest = vectorsPerclass(c(ceil(perclassNum * sampPro) + 1:perclassNum), :);  % select test samples
    testIndexK = index(c(ceil(perclassNum * sampPro) + 1:perclassNum));
    testIndex = [testIndex; testIndexK];
    trainVectors = [trainVectors; selectTrain];
    trainLabels = [trainLabels; repmat(k, ceil(perclassNum * sampPro), 1)];
    testVectors = [testVectors; selectTest];
    testLabels = [testLabels; repmat(k, perclassNum - ceil(perclassNum * sampPro), 1)];
end
[trainVectors, M, m] = scaleData(trainVectors);
[testVectors] = scaleData(testVectors, M, m);

